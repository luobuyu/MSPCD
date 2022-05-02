#include "Config.h"
#include "Instance.h"
#include "Solve.h"
bool SHOW_SOL_LOG = true;
void run_single_instance(Environment& env) {
	//freopen("./log.txt", "w", stderr);
	Solve solve(env, cfg);
	solve.run();
	if (SHOW_SOL_LOG)
	{
		solve.record_sol(env.solution_path_with_time());
		solve.check();
		solve.record_log(env.log_path());
	}
	else
	{
		solve.record_sol();
		solve.record_log();
	}
}


int main(int argc, char* argv[]) {
	Environment env;
	if (argc == 3) // exe  timeout  seed  machine_num
	{
		cfg.timeout = atoi(argv[1]);
		cfg.random_seed = atoi(argv[2]);
		SHOW_SOL_LOG = false;
	}
	else if (argc == 4) // exe  timeout  seed  instance_name
	{
		// instance_name 格式 orga<n>_<r>_<p> 
		// 读取出来 p 
		cfg.timeout = atoi(argv[1]);
		cfg.random_seed = atoi(argv[2]);
		env.set_instance_name(argv[3], cfg.machine_num);
	}
	else 
	{
		cfg.timeout = 60 * 10;
		//cfg.timeout = 10;
		cfg.random_seed = 100;
		string in_name = "ogra500_90_4"; //
		env.set_instance_name(in_name, cfg.machine_num);
		freopen(env.instance_path().c_str(), "r", stdin);
	}
	run_single_instance(env);
	return 0;
}
