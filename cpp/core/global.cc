#include "global.h"
#include <iostream>
using namespace std;
int _my_rank = 0;
int _num_workers = 0;

DEFINE_string(vfile, "../data/test.v", "Location of vertex file");
DEFINE_string(efile, "../data/test.e", "Location of edge file");
DEFINE_string(query_dir, "../data/query/", "Location of query files");
DEFINE_string(result_dir, "../result/test.result", "Dir of result file");
DEFINE_int32(file_location, 0,
                "File location: 0-local(default), 1-HDFS, 2-AWS S3, 5-HTTP");
DEFINE_int32(query_index, 1, "query_index");
DEFINE_string(view_file, "", "Location of view file");
DEFINE_string(r_file, "", "Location of partition file");
DEFINE_string(base_qfile, "", "Location of base query file");
DEFINE_string(base_add_file, "", "Location of base_add_file");
DEFINE_string(base_remove_file, "", "Location of base_remove_file");

