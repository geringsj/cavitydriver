/**
 * taken from 
 * http://c.learncodethehardway.org/book/ex20.html
 *
 * those defines are very useful, but we did not invent them!
 */

#ifndef __dbg_h__
#define __dbg_h__

#include <stdio.h>
#include <errno.h>
#include <string.h>

#ifdef NDEBUG
	#define debug(M, ...)
#else
	#define debug(M, ...) fprintf(stdout/*stderr*/, "[DEBUG] " M " (%s:%d)\n", ##__VA_ARGS__, __FILE__, __LINE__);
#endif

#define clean_errno() (errno == 0 ? "None" : strerror(errno))

#define log_err(M, ...) fprintf(stdout/*stderr*/, "[ERROR] " M " (%s:%d:errno:%s)\n", ##__VA_ARGS__, __FILE__, __LINE__, clean_errno())

#define log_warn(M, ...) fprintf(stdout/*stderr*/, "[WARN] " M " (%s:%d:errno:%s)\n", ##__VA_ARGS__, __FILE__, __LINE__, clean_errno())

#define log_info(M, ...) fprintf(stdout, "[INFO] " M "\n", ##__VA_ARGS__)

#define check(A, M, ...) if(!(A)) { log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define sentinel(M, ...)  { log_err(M, ##__VA_ARGS__); errno=0; goto error; }

#define check_mem(A) check((A), "Out of memory.")

#define check_debug(A, M, ...) if(!(A)) { debug(M, ##__VA_ARGS__); errno=0; goto error; }

#endif
