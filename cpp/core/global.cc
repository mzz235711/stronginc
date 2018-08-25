#include "global.h"
#include <iostream>
using namespace std;
int _my_rank = 0;
int _num_workers = 0;

//DEFINE_string(vfile, "../data/test.v", "Location of vertex file");
//DEFINE_string(efile, "../data/test.e", "Location of edge file");
//DEFINE_string(query_dir, "../data/query/", "Location of query files");
//DEFINE_string(result_dir, "../result/test.result", "Dir of result file");
//DEFINE_int32(file_location, 0,
//                "File location: 0-local(default), 1-HDFS, 2-AWS S3, 5-HTTP");
DEFINE_string(vfile, "../data/test.v", "Location of vertex file");
DEFINE_string(efile, "../data/test.e", "Location of edge file");
DEFINE_string(viewfile, "../data/test.view", "Location of view file");
DEFINE_string(rfile, "../data/test.r", "Location of partition file");
DEFINE_string(query_file, "../data/query/query", "Location of query file");
DEFINE_string(base_add_file, "../data/inc/add_e", "Location of add edge file");
DEFINE_string(base_remove_file, "../data/inc/rm_e", "Location of remove edge file");
DEFINE_int32(query_num, 1, "Query number");
DEFINE_int32(extend_num, 0, "extend file number");
