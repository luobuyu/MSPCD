#ifndef _SRC_CONFIG_H_
#define _SRC_CONFIG_H_
#include <random>
static constexpr int INF = 0x3f3f3f3f;
static constexpr long long INF_LL = 0x3f3f3f3f3f3f3f3f;
static constexpr long long mod = 1e9 + 7;
static constexpr double eps = 1e-5;


using coord_t = int;


static const char* ins_list[]{
	"",
	"",
	""
};

static const int D[][8] = {
	{0, 1, 1, 2, 1, 2, 2, 3},
	{1, 0, 2, 1, 2, 1, 3, 2},
	{1, 2, 0, 1, 2, 3, 1, 2},
	{2, 1, 1, 0, 3, 2, 2, 1},
	{1, 2, 2, 3, 0, 1, 1, 2},
	{2, 1, 3, 2, 1, 0, 2, 1},
	{2, 3, 1, 2, 1, 2, 0, 1},
	{3, 2, 2, 1, 2, 1, 1, 0}
};

struct Config {
	unsigned int random_seed = 100;
	int timeout = 3600;  // ³¬Ê±Ê±¼ä
	int machine_num = 2;
	int ccr = 1;
	int max_tabu_times = 1000;
	int randmove_max_iter = 10;
	Config():timeout(3600)
	{
		random_seed = std::random_device{}();
	}
	Config(int timeout): timeout(timeout)
	{
		random_seed = std::random_device{}();
	}
	Config(int timeout, unsigned int seed) : timeout(timeout), random_seed(seed) {};
} cfg;
#endif // !_SRC_CONFIG_H_