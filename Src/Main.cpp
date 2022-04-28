#include "Config.h"
#include "Instance.h"
#include "Solve.h"
bool SHOW_SOL_LOG = true;
void run_single_instance(Environment& env) {
	Solve solve(env, cfg);
	solve.run();
	if (SHOW_SOL_LOG)
	{
		//solve.record_sol(env.solution_path_with_time());
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
	if (argc == 3) // exe  timeout  seed
	{
		cfg.timeout = atoi(argv[1]);
		cfg.random_seed = atoi(argv[2]);
		SHOW_SOL_LOG = false;
	}
	else if (argc == 4) // exe  timeout  seed  instance_name
	{
		cfg.timeout = atoi(argv[1]);
		cfg.random_seed = atoi(argv[2]);
		env.set_instance_name(argv[3]);
	}
	else 
	{
		cfg.timeout = 60;
		cfg.random_seed = 10;
		string in_name = "t50_40_6";
		env.set_instance_name(in_name);
		freopen(env.instance_path().c_str(), "r", stdin);
	}
	run_single_instance(env);
	return 0;
}
