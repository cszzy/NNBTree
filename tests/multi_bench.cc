#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <thread>
#include <getopt.h>
#include <unistd.h>
#include <cmath>
#include <vector>
#include <set>
#include <map>

#include "db_interface.h"
#include "timer.h"
#include "util.h"
#include "random.h"

#include "apex_util.h"

using ycsbc::KvDB;

bool static_lru;
uint64_t miss_times[64];
uint64_t evict_times[64];
// std::unordered_set<char*> subtree_set[64];

using namespace util;

// 实时获取程序占用的内存，单位：kb
size_t physical_memory_used_by_process()
{
    FILE* file = fopen("/proc/self/status", "r");
    int result = -1;
    char line[128];

    while (fgets(line, 128, file) != nullptr) {
        if (strncmp(line, "VmRSS:", 6) == 0) {
            int len = strlen(line);

            const char* p = line;
            for (; std::isdigit(*p) == false; ++p) {}

            line[len - 3] = 0;
            result = atoi(p);

            break;
        }
    }

    fclose(file);
    return result;
}

void clear_cache() {
  // Remove cache
  int size = 256 * 1024 * 1024;
  char *garbage = new char[size];
  for (int i = 0; i < size; ++i)
    garbage[i] = i;
  for (int i = 100; i < size; ++i)
    garbage[i] += garbage[i - 100];
  delete[] garbage;
}

struct operation
{
  /* data */
  uint32_t op_type;
  uint32_t op_len; // for only scan
  uint64_t key_num;
};

void show_help(char* prog) {
  std::cout <<
    "Usage: " << prog << " [options]" << std::endl <<
    std::endl <<
    "  Option:" << std::endl <<
    "    --thread[-t]             thread number" << std::endl <<
    "    --load-size              LOAD_SIZE" << std::endl <<
    "    --put-size               PUT_SIZE" << std::endl <<
    "    --get-size               GET_SIZE" << std::endl <<
    "    --workload               WorkLoad" << std::endl <<
    "    --help[-h]               show help" << std::endl;
}

int thread_num = 32;
size_t LOAD_SIZE   = 150000000;
size_t PUT_SIZE    = 6000000;
size_t GET_SIZE    = 1000000;
size_t DELETE_SIZE = 1000000;
int Loads_type = 0;


uint64_t get_cellid(const std::string &line) {
  uint64_t id;
  double lat, lon;
  std::stringstream strin(line);
          strin >> id >> lon >> lat;
  return id;
}

double get_longitude(const std::string &line) {
  uint64_t id;
  double lat, lon;
  std::stringstream strin(line);
          strin >> id >> lon >> lat;
  return lon;
}

double get_lat(const std::string &line) {
  uint64_t id;
  double lat, lon;
  std::stringstream strin(line);
          strin >> id >> lon >> lat;
  return lat;
}

uint64_t get_longlat(const std::string &line) {
  uint64_t id;
  double lat, lon;
  std::stringstream strin(line);
          strin >> id >> lon >> lat;
  return (lon * 180 + lat) * 1e7;
}

template<typename T>
std::vector<T>read_data_from_osm(const std::string load_file, 
    T (*get_data)(const std::string &) = []{ return static_cast<T>(0);},
    const std::string output = "/home/sbh/generate_random_osm_longlat.dat")
{
  std::vector<T> data;
  std::set<T> unique_keys;
  std::cout << "Use: " << __FUNCTION__ << std::endl;
    const uint64_t ns = util::timing([&] { 
      std::ifstream in(load_file);
      if (!in.is_open()) {
        std::cerr << "unable to open " << load_file << std::endl;
        exit(EXIT_FAILURE);
      }
      uint64_t size = 0;
      while (!in.eof())
      {
        /* code */
        std::string tmp;
        getline(in, tmp); // 去除第一行
        while(getline(in, tmp)) {
          T key = get_data(tmp);
          unique_keys.insert(key);
          size ++;
          if(size % 100000 == 0) std::cerr << "Load: " << size << "\r";
        }
      }
      in.close();
      std::cerr << "Finshed loads ......\n";
      data.assign(unique_keys.begin(), unique_keys.end());
      std::random_shuffle(data.begin(), data.end());
      size = data.size();
      std::cerr << "Finshed random ......\n"; 
      std::ofstream out(output, std::ios::binary);
      out.write(reinterpret_cast<char*>(&size), sizeof(uint64_t));
      out.write(reinterpret_cast<char*>(data.data()), data.size() * sizeof(uint64_t));
      out.close(); 
      std::cout << "read size: " << size << ", unique data: " << unique_keys.size() << std::endl;
  });
  const uint64_t ms = ns/1e6;
  std::cout << "generate " << data.size() << " values in "
            << ms << " ms (" << static_cast<double>(data.size())/1000/ms
            << " M values/s)" << std::endl;   
  return data;
}

template<typename T>
std::vector<T>load_data_from_osm(
    const std::string dataname = "/home/sbh/generate_random_osm_cellid.dat")
{
  return util::load_data<T>(dataname);
}

std::vector<uint64_t> generate_random_ycsb(size_t op_num)
{
  std::vector<uint64_t> data; 
  data.resize(op_num);
  std::cout << "Use: " << __FUNCTION__ << std::endl;
  const uint64_t ns = util::timing([&] { 
    nnbtree::Random rnd(0, op_num - 1);
    for (size_t i = 0; i < op_num; ++i)
      data[i] = utils::Hash(i);
    // for (size_t i = 0; i < op_num; ++i)
    //   std::swap(data[i], data[rnd.Next()]);
  });
  const uint64_t ms = ns/1e6;
  std::cout << "generate " << data.size() << " values in "
            << ms << " ms (" << static_cast<double>(data.size())/1000/ms
            << " M values/s)" << std::endl;   
  return data;
}

std::vector<uint64_t> generate_uniform_random(size_t op_num)
{
  std::vector<uint64_t> data; 
  data.resize(op_num);
  std::cout << "Use: " << __FUNCTION__ << std::endl;
  const uint64_t ns = util::timing([&] { 
    nnbtree::Random rnd(0, UINT64_MAX);
    for (size_t i = 0; i < op_num; ++i)
      data[i] = rnd.Next();
  });
  const uint64_t ms = ns/1e6;
  std::cout << "generate " << data.size() << " values in "
            << ms << " ms (" << static_cast<double>(data.size())/1000/ms
            << " M values/s)" << std::endl;  
  return data;
}

int main(int argc, char *argv[]) {
    static struct option opts[] = {
  /* NAME               HAS_ARG            FLAG  SHORTNAME*/
    {"thread",          required_argument, NULL, 't'},
    {"load-size",       required_argument, NULL, 0},
    {"put-size",        required_argument, NULL, 0},
    {"get-size",        required_argument, NULL, 0},
    {"dbname",          required_argument, NULL, 0},
    {"workload",        required_argument, NULL, 0},
    {"loadstype",       required_argument, NULL, 0},
    {"help",            no_argument,       NULL, 'h'},
    {NULL, 0, NULL, 0}
  };

  int c;
  int opt_idx;
  std::string  dbName= "nnbtree";
  std::string  load_file= "";
  while ((c = getopt_long(argc, argv, "t:s:dh", opts, &opt_idx)) != -1) {
    switch (c) {
      case 0:
        switch (opt_idx) {
          case 0: thread_num = atoi(optarg); break;
          case 1: LOAD_SIZE = atoi(optarg); break;
          case 2: PUT_SIZE = atoi(optarg); break;
          case 3: GET_SIZE = atoi(optarg); break;
          case 4: dbName = optarg; break;
          case 5: load_file = optarg; break;
          case 6: Loads_type = atoi(optarg); break;
          case 7: show_help(argv[0]); return 0;
          default: std::cerr << "Parse Argument Error!" << std::endl; abort();
        }
        break;
      case 't': thread_num = atoi(optarg); break;
      case 'h': show_help(argv[0]); return 0;
      case '?': break;
      default:  std::cout << (char)c << std::endl; abort();
    }
  }

  std::cout << "THREAD NUMBER:         " << thread_num << std::endl;
  std::cout << "LOAD_SIZE:             " << LOAD_SIZE << std::endl;
  std::cout << "PUT_SIZE:              " << PUT_SIZE << std::endl;
  std::cout << "GET_SIZE:              " << GET_SIZE << std::endl;
  std::cout << "DB  name:              " << dbName << std::endl;
  std::cout << "Workload:              " << load_file << std::endl;

  std::vector<uint64_t> data_base;
  switch (Loads_type)
  {
  case -2:
    data_base = read_data_from_osm<uint64_t>(load_file, get_cellid);
    break;
  case -1:
    data_base = read_data_from_osm<uint64_t>(load_file, get_longlat);
    break;
  case 0:
    data_base = generate_uniform_random(LOAD_SIZE + PUT_SIZE * 8);
    break;
  case 1:
    data_base = generate_random_ycsb(LOAD_SIZE + PUT_SIZE * 8);
    break;
  case 2:
    data_base = load_data_from_osm<uint64_t>("/home/sbh/generate_random_osm_cellid.dat");;
    break;
  case 3:
    data_base = load_data_from_osm<uint64_t>("/home/zzy/dataset/generate_random_ycsb.dat");;
    break;
  case 4:
    data_base = load_data_from_osm<uint64_t>("/home/zzy/dataset/generate_random_osm_longlat.dat");;
    break;
  case 5:
    data_base = load_data_from_osm<uint64_t>("/home/zzy/dataset/generate_random_osm_longtitudes.dat");
    break;
  case 6:
    data_base = load_data_from_osm<uint64_t>("/home/zzy/dataset/lognormal.dat");
    break;
  default:
    data_base = generate_uniform_random(LOAD_SIZE + PUT_SIZE * 8);
    break;
  }

  // NVM::env_init();
  nnbtree::init_numa_pool();
  
  KvDB *db = nullptr;
  if (dbName == "fastfair")
  {
    db = new nnbtree::fastfairDB();
  } else if (dbName == "nnbtree") {
    db = new nnbtree::nnbtreeDB();
  } else {
    assert(false);
  }

  db->Init();

  GET_SIZE = 64000000;
  uint64_t *GET_data = apex::get_search_keys_zipf_with_theta<uint64_t>(data_base.data(), LOAD_SIZE + PUT_SIZE, GET_SIZE, 0.99);
  // uint64_t *GET_data = apex::get_search_keys<uint64_t>(data_base.data(), LOAD_SIZE + PUT_SIZE, GET_SIZE);

  // {
  //   unordered_set<uint64_t> sss;
  //   for (int i = 0; i < GET_SIZE; i++) {
  //     sss.insert(GET_data[i]);
  //   }

  //   std::cout << "zipfan key nums: " << sss.size() << std::endl;
  // }

  size_t init_dram_space_use = physical_memory_used_by_process();
  std::cout << "before newdb, dram space use: " << init_dram_space_use / 1024.0 /1024.0  << " GB" << std::endl;

  nnbtree::Timer timer;
  uint64_t us_times; 
  std::cout << "Start run ...." << std::endl;
  // {
  //   int init_size = LOAD_SIZE;
  //   std::mt19937_64 gen_payload(std::random_device{}());
  //   auto values = new std::pair<uint64_t, uint64_t>[init_size];
  //   for (int i = 0; i < init_size; i++) {
  //     // values[i].first = data_base[data_base.size() - i - 1];
  //     values[i].first = data_base[i];
  //     values[i].second = static_cast<uint64_t>(gen_payload());
  //   }
  //   std::sort(values, values + init_size,
  //             [](auto const& a, auto const& b) { return a.first < b.first; });
  //   db->Bulk_load(values, init_size);
  // }
  static_lru = false;
  memset(miss_times, 0, sizeof(miss_times));
  memset(evict_times, 0, sizeof(evict_times));
  {
     // Load
    timer.Record("start");
    std::vector<std::thread> threads;
    std::atomic_int thread_id_count(0);
    size_t per_thread_size = LOAD_SIZE / thread_num;
    for(int i = 0; i < thread_num; i ++) {
        threads.emplace_back([&](){
            int thread_id = thread_id_count.fetch_add(1);
            nnbtree::my_thread_id = thread_id;
            // std::cout << thread_id << std::endl << std::flush;
#ifdef NUMA_TEST
            nnbtree::bindCore(nnbtree::my_thread_id);
#endif
            size_t start_pos = thread_id * per_thread_size;
            size_t size = (thread_id == thread_num-1) ? LOAD_SIZE-(thread_num-1)*per_thread_size : per_thread_size;
            for (size_t j = 0; j < size; ++j) {
                auto ret = db->Put(data_base[start_pos+j], data_base[start_pos+j]);
                if (ret != 1) {
                    std::cout << "load error, key: " << data_base[start_pos+j] << ", size: " << j << std::endl;
                    db->Put(data_base[start_pos+j], data_base[start_pos+j]);
                    assert(0);
                }

                if(thread_id == 0 && (j + 1) % 1000000 == 0) {
                  std::cout << "Operate: " << j + 1 << '\r' << std::endl;  
                  std::cout << "dram space use: " << (physical_memory_used_by_process() - init_dram_space_use) / 1024.0 /1024.0  << " GB" << std::endl;
                  db->Info();
                } 
            }
        });
    }
    for (auto& t : threads) {
      t.join();
    }

    timer.Record("stop");
    us_times = timer.Microsecond("stop", "start");
    std::cout << "[Metic-Load]: Load " << LOAD_SIZE << ": " 
              << "cost " << us_times/1000000.0 << "s, " 
              << "iops " << (double)(LOAD_SIZE)/(double)us_times*1000000.0 << " ." << std::endl;
  }

  static_lru = false;
  memset(miss_times, 0, sizeof(miss_times));
  memset(evict_times, 0, sizeof(evict_times));
  {
     // Put
    clear_cache();
    std::vector<std::thread> threads;
    std::atomic_int thread_id_count(0);
    size_t per_thread_size = PUT_SIZE / thread_num;
    timer.Clear();
    timer.Record("start");
    for(int i = 0; i < thread_num; i ++) {
        threads.emplace_back([&](){
            int thread_id = thread_id_count.fetch_add(1);
            nnbtree::my_thread_id = thread_id;
#ifdef NUMA_TEST
            nnbtree::bindCore(nnbtree::my_thread_id);
#endif
            size_t start_pos = thread_id * per_thread_size + LOAD_SIZE;
            size_t size = (thread_id == thread_num-1) ? PUT_SIZE-(thread_num-1)*per_thread_size : per_thread_size;
            for (size_t j = 0; j < size; ++j) {
                auto ret = db->Put(data_base[start_pos+j], data_base[start_pos+j]);
                if (ret != 1) {
                    std::cout << "Put error, key: " << data_base[start_pos+j] << ", size: " << j << std::endl;
                    assert(0);
                }
                if(thread_id == 0 && (j + 1) % 100000 == 0) std::cerr << "Operate: " << j + 1 << '\r'; 
            }
        });
    }
    for (auto& t : threads) {
      t.join();
    }
        
    timer.Record("stop");
    us_times = timer.Microsecond("stop", "start");
    std::cout << "[Metic-Put]: Put " << PUT_SIZE << ": " 
              << "cost " << us_times/1000000.0 << "s, " 
              << "iops " << (double)(PUT_SIZE)/(double)us_times*1000000.0 << " ." << std::endl;
  }
  // std::cout << "getchar:" <<std::endl;
  // getchar();
  static_lru = true;
  memset(miss_times, 0, sizeof(miss_times));
  memset(evict_times, 0, sizeof(evict_times));
  // GET_SIZE = 100000000;
  {
     // Get
    clear_cache();
    // std::cout << "getchar" << std::endl;
    // getchar();

    std::vector<std::thread> threads;
    std::atomic_int thread_id_count(0);
    size_t per_thread_size = GET_SIZE / thread_num;

    timer.Clear();
    timer.Record("start");
    
    for (int i = 0; i < thread_num; ++i) {
        threads.emplace_back([&](){
            thread_local uint64_t error_gets = 0;
            int thread_id = thread_id_count.fetch_add(1);
            nnbtree::my_thread_id = thread_id;
#ifdef NUMA_TEST
            nnbtree::bindCore(nnbtree::my_thread_id);
#endif
            size_t start_pos = thread_id *per_thread_size;
            size_t size = (thread_id == thread_num-1) ? GET_SIZE-(thread_num-1)*per_thread_size : per_thread_size;
            size_t value;
            for (size_t j = 0; j < size; ++j) {
                bool ret = db->Get(GET_data[start_pos+j], value);
                if (ret != true || value != GET_data[start_pos+j]) {
                    // std::cout << "Get error!" << std::endl;
                    error_gets++;
                }
                if(thread_id == 0 && (j + 1) % 100000 == 0) std::cerr << "Operate: " << j + 1 << '\r'; 
            }
            if (thread_id == 0)
              std::cout << "eror_gets: " << error_gets << std::endl;
        });
    }

    for (auto& t : threads) {
      t.join();
    }

    timer.Record("stop");
    us_times = timer.Microsecond("stop", "start");

    uint64_t total_miss_times = 0;
    uint64_t total_evict_times = 0;
    for (int i = 0; i < 64; i++) {
      total_miss_times += miss_times[i];
      total_evict_times += evict_times[i];
    }

    // std::unordered_set<char*> sssss;
    // for (int i = 0; i < 64; i++) {
    //   for(auto iter = subtree_set->begin(); iter != subtree_set->end(); iter++) {
    //     sssss.insert(*iter);
    //   }
    // }
    
    std::cout << "total_miss_times:" << total_miss_times << ", total_evict_times: " << total_evict_times << std::endl;
    // std::cout << "subtree set size: " << sssss.size() << std::endl;
    std::cout << "[Metic-Get]: Get " << GET_SIZE << ": " 
              << "cost " << us_times/1000000.0 << "s, " 
              << "iops " << (double)(GET_SIZE)/(double)us_times*1000000.0 << " ." << std::endl;
  }

  // mixed test
  GET_SIZE = 10000000;
  {
    // std::vector<float> insert_ratios = {0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    std::vector<float> insert_ratios = {0.3, 0.7};
    float insert_ratio = 0;
    for (int ini = 0; ini < insert_ratios.size(); ini++) {
      insert_ratio = insert_ratios[ini];
      // uint64_t *mix_test = apex::get_search_keys_zipf_with_theta<uint64_t>(data_base.data(), LOAD_SIZE + PUT_SIZE, GET_SIZE, 0.9);
      uint64_t *mix_test = apex::get_search_keys<uint64_t>(data_base.data(), LOAD_SIZE + PUT_SIZE, GET_SIZE);

      std::vector<std::thread> threads;
      std::atomic_int thread_id_count(0);
      size_t per_thread_size = GET_SIZE / thread_num;

      vector<vector<float>> insert_ratio_arr(thread_num, vector<float>());
      util::FastRandom ranny(18);
      for (int i = 0; i < thread_num; i++) {
        for (int j = 0; j < GET_SIZE / (thread_num == 1 ? 1 : thread_num - 1); j++) {
          insert_ratio_arr[i].push_back(ranny.ScaleFactor());
        }
      }

      clear_cache();

      timer.Clear();
      timer.Record("start");
      for (int i = 0; i < thread_num; ++i) {
          threads.emplace_back([&](){
              int thread_id = thread_id_count.fetch_add(1);
              size_t get_start_pos = thread_id *per_thread_size;
              size_t put_start_pos = thread_id * per_thread_size + LOAD_SIZE + PUT_SIZE * (1 + ini);
              size_t size = (thread_id == thread_num-1) ? GET_SIZE-(thread_num-1)*per_thread_size : per_thread_size;
              size_t value;
              int put_op = 0, get_op = 0;
              for (size_t j = 0; j < size; ++j) {
                  if (insert_ratio_arr[thread_id][j] < insert_ratio) {
                    put_op++;
                    int ret = db->Put(data_base[put_start_pos], data_base[put_start_pos++]);
                    if (ret != 1) {
                      std::cout << "Put error, key: " << data_base[put_start_pos-1] << ", size: " << j << std::endl;
                      assert(0);
                    }
                  } else {
                    get_op++;
                    int ret = db->Get(mix_test[get_start_pos+j], value);
                    if (ret != 1 || value != mix_test[get_start_pos+j]) {
                        // std::cout << "Get error!" << std::endl;
                    }
                  }
                  if(thread_id == 0 && (j + 1) % 100000 == 0) {
                    std::cout << "Operate: " << j + 1 << " get_op: " << get_op << " put_op: " << put_op << std::endl;
                    put_op = get_op = 0;
                  }  
              }
          });
      }

      for (auto& t : threads) {
        t.join();
      }
          
      timer.Record("stop");
      us_times = timer.Microsecond("stop", "start");
      std::cout << "[Metic-Mix]: Mix "<< insert_ratio << " " << GET_SIZE
                << "cost " << us_times/1000000.0 << "s, " 
                << "iops " << (double)(GET_SIZE)/(double)us_times*1000000.0 << " ." << std::endl;
    }
  }

  delete db;

  NVM::env_exit();
  return 0;
}