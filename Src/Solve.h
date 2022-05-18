#ifndef _SRC_SOLVE_H_
#define _SRC_SOLVE_H_

#include "Instance.h"
#include "Utils.h"
#include <queue>
#include <set>
#include <cassert>
using namespace utils;
using namespace fast_io;
using namespace Debug;
//#define SHOW_MOVE
class Solve {

private:
	const Environment& _env;
	const Config& _cfg;

	Instance _ins;
	default_random_engine _gen;
	Timer timer;
	vector<vector<Operation>> machine;		// ������ݽṹ
	vector<vector<Operation>> best_machine;	// ���Ž�	
	vector<pair<int, int> > task_of_machine;
	vector<int> R;	// start �� i ���·
	vector<int> Q;	// end �� i ���·
	vector<vector<Block>> critical_block;
	
	// �������rowһ��Ԫ��
	vector<MoveSeq> best_tabu_move_seqs;
	vector<MoveSeq> best_not_tabu_move_seqs;
	vector<MoveMachine> best_tabu_move_machines;
	vector<MoveMachine> best_not_tabu_move_machines;
	vector<MoveMachine> best_move_machine_score_eqs;
	// ���ɱ�ṹ
	TabuTable tabu_list;
	TabuArr tabu_machine; // task_id, machine_id;

	int makespan;	  // �Ż�Ŀ��
	int best_makespan;// ���Ž�

	int Start, End;

	double _duration; // ���Ž����ʱ��
	int _iteration;   // ���Ž���ֵ�������
	int _all_iteration; // �ܹ��ĵ�������
public:

	Solve() = delete;
	
	Solve(const Environment& env, Config& cfg) :
		_env(env), _cfg(cfg), _ins(env), _gen(_cfg.random_seed),
		timer(cfg.timeout), _duration(0), _iteration(0),
		R(_ins.get_vir_task_num()), Q(_ins.get_vir_task_num()),
		Start(0), End(_ins.get_task_num() + 1), task_of_machine(_ins.get_vir_task_num()),
		makespan(INF), best_makespan(INF),
		machine(_ins.get_machine_num()), best_machine(_ins.get_machine_num()),
		critical_block(_ins.get_machine_num()), tabu_list(_ins.get_machine_num(), 5, 10),
		tabu_machine(_ins.get_vir_task_num(), _ins.get_machine_num(), 20, 30),
		_all_iteration(0)
		 {
	}	// �� _obj_area��ʼ��Ϊ���ֵ

	// [lb, ub)
	int fast_rand(int lb, int ub)
	{
		// [0, ub - lb - 1] -> [lb, ub - 1]
		return _gen() % (ub - lb) + lb;
	}
	// [0, ub)
	int fast_rand(int ub)
	{
		// [0, ub - 1]
		return _gen() % ub;
	}

	void sumR()
	{
		queue<int> q;
		vector<int> in_deg(_ins.get_vir_task_num(), 0);
		
		// ��һ�����
		for (int i = 0; i < _ins.get_vir_task_num(); ++i) { in_deg[i] = _ins.get_jp()[i].size(); }
		for (int i = 0; i < machine.size(); ++i)
		{
			for (int j = 0; j < machine[i].size(); ++j)
			{
				in_deg[machine[i][j].task_id] += (j != 0); // ֻҪ���ǻ��������еĵ�һ��Ԫ�ؾ���Ҫ��mp
				if (in_deg[machine[i][j].task_id] == 0) machine[i][j].R = 0;
				else machine[i][j].R = -1;
			}
		}

		q.push(Start);
		while (!q.empty())
		{
			int u = q.front(); q.pop();
			auto& task_u = _ins.get_tasks()[u];
			int mp = machine[task_u.belong_machine][task_u.machine_seq].get_r_time();
			// ���һ��������еĺ���
			if (task_u.machine_seq + 1 < machine[task_u.belong_machine].size())
			{
				auto& u_next_op = machine[task_u.belong_machine][task_u.machine_seq + 1];
				u_next_op.R = max(mp, u_next_op.R);
				//in_deg[machine[task_u.belong_machine][task_u.machine_seq + 1].task_id]--;
				if (--in_deg[u_next_op.task_id] == 0)
				{
					q.push(u_next_op.task_id);
				}
			}
			// ����task�����еĺ���
			for (int i = 0; i < _ins.get_js()[u].size(); ++i)
			{
				auto& edge = _ins.get_js()[u][i];
				auto& js_u_task = _ins.get_tasks()[edge.v];
				int cc_time = D[task_u.belong_machine][js_u_task.belong_machine] * edge.w * cfg.ccr;
				machine[js_u_task.belong_machine][js_u_task.machine_seq].R = max(machine[js_u_task.belong_machine][js_u_task.machine_seq].R, mp + cc_time);
				if (--in_deg[edge.v] == 0)
				{
					q.push(edge.v);
				}
			}
		}
	}

	void sumQ()
	{
		queue<int> q;
		vector<int> in_deg(_ins.get_vir_task_num(), 0);

		// ��һ����ȣ���ʵ�ǳ���
		for (int i = 0; i < _ins.get_vir_task_num(); ++i) { in_deg[i] = _ins.get_js()[i].size(); }
		for (int i = 0; i < machine.size(); ++i)
		{
			int len = machine[i].size();
			for (int j = 0; j < len; ++j)
			{
				in_deg[machine[i][j].task_id] += (j != len - 1);
				if (in_deg[machine[i][j].task_id] == 0) machine[i][j].Q = 0;
				else machine[i][j].Q = -1;
			}
		}

		q.push(End);
		while (!q.empty())
		{
			int u = q.front(); q.pop();
			//cerr << u << " ";
			auto& task_u = _ins.get_tasks()[u];
			int ms = machine[task_u.belong_machine][task_u.machine_seq].get_q_time();
			// ���һ��������еĺ���
			if (task_u.machine_seq - 1 >= 0)
			{
				auto& u_pre_op = machine[task_u.belong_machine][task_u.machine_seq - 1];
				u_pre_op.Q = max(ms, u_pre_op.Q);
				//in_deg[machine[task_u.belong_machine][task_u.machine_seq + 1].task_id]--;
				if (--in_deg[u_pre_op.task_id] == 0)
				{
					q.push(u_pre_op.task_id);
				}
			}
			// ����task�����е�ǰ׺
			for (int i = 0; i < _ins.get_jp()[u].size(); ++i)
			{
				auto& edge = _ins.get_jp()[u][i];
				auto& jp_u_task = _ins.get_tasks()[edge.v];
				int cc_time = D[task_u.belong_machine][jp_u_task.belong_machine] * edge.w * cfg.ccr;
				machine[jp_u_task.belong_machine][jp_u_task.machine_seq].Q = max(machine[jp_u_task.belong_machine][jp_u_task.machine_seq].Q, ms + cc_time);
				if (--in_deg[edge.v] == 0)
				{
					q.push(edge.v);
				}
			}
		}
		//cerr << endl;
	}

	void sum_makespan()
	{
		auto& start_task = _ins.get_tasks()[Start];
		auto& start_op = machine[start_task.belong_machine][start_task.machine_seq];
		auto& end_task = _ins.get_tasks()[End];
		auto& end_op = machine[end_task.belong_machine][end_task.machine_seq];
		assert(start_op.get_r_q_time() == end_op.get_r_q_time());
		makespan = start_op.get_r_q_time();
		assert(makespan != -1);
	}

	void sum_critical_path()
	{
		sumR();
		sumQ();
		sum_makespan();
		for (int i = 0; i < critical_block.size(); ++i) { critical_block[i].clear(); }
		for (int i = 0; i < machine.size(); ++i)
		{
			for (int j = 0; j < machine[i].size();)
			{
				if (machine[i][j].get_r_q_time() == makespan 
					&& machine[i][j].task_id != Start 
					&& machine[i][j].task_id != End)
				{
					int start_index = j;
					int end_index = j;
					j++;
					while (j < machine[i].size() && machine[i][j].get_r_q_time() == makespan
						&& machine[i][j].task_id != Start && machine[i][j].task_id != End)
					{
						end_index = j;
						j++;
					}
					critical_block[i].push_back({ start_index, end_index });
				}
				else
				{
					j++;
				}
			}
		}
	}

	void init_solution()
	{
		// Ϊÿһ�������������һ������
		// ��������һ��������ȽϺ�
		queue<int> q;
		vector<int> in_deg(_ins.get_vir_task_num(), 0);
		for (int i = 0; i < _ins.get_vir_task_num(); ++i) { in_deg[i] = _ins.get_jp()[i].size(); }
		q.push(Start);
		while (!q.empty())
		{
			int u = q.front(); q.pop();
			/**************�Ӷ�����ȡ��������һ̨����*****************/
			int machine_id = fast_rand(machine.size());
			_ins.get_tasks()[u].belong_machine = machine_id;
			//machine[machine_id].push_back(_ins.get_tasks()[u]);
			for (int i = 0; i < _ins.get_js()[u].size(); ++i)
			{
				int v = _ins.get_js()[u][i].v;
				// ���v�����Ϊ 0����ô�������
				if (--in_deg[v] == 0)
				{
					q.push(v);
				}
			}
		}

		for (int i = 0; i < machine.size(); ++i) machine[i].clear();

		for (int i = 0; i < _ins.get_vir_task_num(); ++i)
		{
			int jp = 0, mp = 0;
			auto& task = _ins.get_tasks()[i];
			int select_machine = task.belong_machine;

			for (int j = 0; j < _ins.get_jp()[i].size(); ++j)
			{
				auto& p_edge = _ins.get_jp()[i][j];
				auto& p_task = _ins.get_tasks()[p_edge.v];
				int cc_time = D[task.belong_machine][p_task.belong_machine] * p_edge.w * cfg.ccr;
				jp = max(jp, p_task.get_total_time() + cc_time);
			}
			
			if (machine[select_machine].size() != 0)
			{
				mp = machine[select_machine].back().R + machine[select_machine].back().dura_time;
			}
			task.start_time = max(jp, mp);
			// operation(int id, int time, int r) :task_id(id), dura_time(time), R(r)
			task.machine_seq = machine[select_machine].size();	// ��Ҫ��¼�ĸ��������ĸ�λ��
			machine[select_machine].push_back(Operation(i, task.length, task.start_time));
		}
	}

	bool check_obj(int iter)
	{
		if (makespan < best_makespan)
		{
			best_makespan = makespan;
			best_machine.clear();
			best_machine = machine;
			for (int i = 0; i < _ins.get_vir_task_num(); ++i)
			{
				task_of_machine[i].first = _ins.get_tasks()[i].belong_machine;
				task_of_machine[i].second = _ins.get_tasks()[i].machine_seq;
			}
			_iteration = iter;	// ��¼���Ž���ֵĵ�������
			_duration = timer.getDuration();	// ��¼����ʱ��
			return true;
		}
		return false;
	}

	// �� u �ƶ��� v ����
	bool is_move_seq_backward_legal(int machine_id, int u, int v)
	{
		auto& task_u_id = machine[machine_id][u].task_id;
		auto& task_v_id = machine[machine_id][v].task_id;
		// ��� js(u) �� End �Ļ���ô��
		for (int i = 0; i < _ins.get_js()[task_u_id].size(); ++i)
		{
			auto& edge = _ins.get_js()[task_u_id][i];
			auto& js_u_task = _ins.get_tasks()[edge.v];
			int cc_time = D[js_u_task.belong_machine][machine_id] * edge.w * cfg.ccr;
			//if ((task_v_id == edge.v) || machine[machine_id][v].Q < machine[js_u_task.belong_machine][js_u_task.machine_seq].Q + cc_time) return false;
			//if ((task_v_id == edge.v) || machine[machine_id][v].get_q_time() < machine[js_u_task.belong_machine][js_u_task.machine_seq].get_q_time() + cc_time) return false;

			if ((task_v_id == edge.v) || machine[machine_id][v].get_q_time() < machine[js_u_task.belong_machine][js_u_task.machine_seq].get_q_time()) return false;
		}
		return true;
	}

	int try_to_move_seq_back(int machine_id, int u, int v)
	{
		// 0, v - u
		int length = v - u + 1;
		int u_task_id = machine[machine_id][u].task_id;
		int v_task_id = machine[machine_id][v].task_id;
		vector<int> R_copy(length), Q_copy(length);
		for (int i = u; i <= v; ++i)
		{
			R_copy[i - u] = machine[machine_id][i].R;
			Q_copy[i - u] = machine[machine_id][i].Q;
		}
		/******************* ����R_copy ************************/
		// R_cpoy[l1]
		int l1 = u + 1; // ��ʾl1�� machine_id �е��±�
		int l1_task_id = machine[machine_id][l1].task_id;
		auto& l1_task = _ins.get_tasks()[l1_task_id];
		int jp = 0, mp = 0;

		for (int i = 0; i < _ins.get_jp()[l1_task_id].size(); ++i)
		{
			auto& edge = _ins.get_jp()[l1_task_id][i];
			auto& jp_l1_task = _ins.get_tasks()[edge.v];
			int cc_time = D[machine_id][jp_l1_task.belong_machine] * edge.w * cfg.ccr;
			jp = max(jp, machine[jp_l1_task.belong_machine][jp_l1_task.machine_seq].get_r_time() + cc_time);
		}
			// ע���Ƿ��ǻ��������еĵ�һ��
		if (u >= 1) mp = machine[machine_id][u - 1].get_r_time();
		R_copy[l1 - u] = max(jp, mp);

		// [l2, v]
		for (int i = u + 2; i <= v; ++i)
		{
			jp = 0;
			// mp(w)
			mp = R_copy[i - 1 - u] + machine[machine_id][i - 1].dura_time;
			// �� jp(w)
			int w_task_id = machine[machine_id][i].task_id;
			for (int j = 0; j < _ins.get_jp()[w_task_id].size(); ++j)
			{
				auto& edge = _ins.get_jp()[w_task_id][j];
				auto& jp_w_task = _ins.get_tasks()[edge.v];
				int cc_time = D[jp_w_task.belong_machine][machine_id] * edge.w * cfg.ccr;
				jp = max(jp, machine[jp_w_task.belong_machine][jp_w_task.machine_seq].get_r_time() + cc_time);
			}
			R_copy[i - u] = max(jp, mp);
		}
		// R_copy[u]
		jp = 0;
		mp = R_copy[v - u] + machine[machine_id][v].dura_time;
		
		for (int i = 0; i < _ins.get_jp()[u_task_id].size(); ++i)
		{
			auto& edge = _ins.get_jp()[u_task_id][i];
			auto& jp_u_task = _ins.get_tasks()[edge.v];
			int cc_time = D[jp_u_task.belong_machine][machine_id] * edge.w * cfg.ccr;
			jp = max(jp, machine[jp_u_task.belong_machine][jp_u_task.machine_seq].get_r_time() + cc_time);
		}
		R_copy[u - u] = max(jp, mp);
		/*************************** ���� Q_copy ********************************/
		// �ȼ��� Q_copy[u]
		int js = 0, ms = 0;
		for (int i = 0; i < _ins.get_js()[u_task_id].size(); ++i)
		{
			auto& edge = _ins.get_js()[u_task_id][i];
			auto& js_u_task = _ins.get_tasks()[edge.v];
			int cc_time = D[js_u_task.belong_machine][machine_id] * edge.w * cfg.ccr;
			js = max(js, machine[js_u_task.belong_machine][js_u_task.machine_seq].get_q_time() + cc_time);
		}
		if (v + 1 <= machine[machine_id].size() - 1)
		{
			ms = machine[machine_id][v + 1].get_q_time();
		}
		Q_copy[u - u] = max(js, ms);
		// �ټ��� Q_copy[v]
		js = 0, ms = 0;
		for (int i = 0; i < _ins.get_js()[v_task_id].size(); ++i)
		{
			auto& edge = _ins.get_js()[v_task_id][i];
			auto& js_v_task = _ins.get_tasks()[edge.v];
			int cc_time = D[js_v_task.belong_machine][machine_id] * edge.w * cfg.ccr;
			js = max(js, machine[js_v_task.belong_machine][js_v_task.machine_seq].get_q_time() + cc_time);
		}
		ms = machine[machine_id][u].dura_time + Q_copy[u - u];
		Q_copy[v - u] = max(js, ms);

		// [l1, lk]
		for (int i = v - 1; i >= u + 1; --i)
		{
			js = 0, ms = 0;
			ms = machine[machine_id][i + 1].dura_time + Q_copy[i - u + 1];
			int w_task_id = machine[machine_id][i].task_id;
			for (int j = 0; j < _ins.get_js()[w_task_id].size(); ++j)
			{
				auto& edge = _ins.get_js()[w_task_id][j];
				auto& js_w_task = _ins.get_tasks()[edge.v];
				int cc_time = D[js_w_task.belong_machine][machine_id] * edge.w * cfg.ccr;
				js = max(js, machine[js_w_task.belong_machine][js_w_task.machine_seq].get_q_time() + cc_time);
			}
			Q_copy[i - u] = max(js, ms);
		}

		// �������֮��Ĵ���
		int new_makespan = 0;
		for (int i = u; i <= v; ++i)
		{
			new_makespan = max(new_makespan, R_copy[i - u] + Q_copy[i - u] + machine[machine_id][i].dura_time);
		}
		return new_makespan;
	}

	// �� v �Ƶ� u ��ǰ��
	bool is_move_seq_forward_legal(int machine_id, int u, int v)
	{
		auto& task_v_id = machine[machine_id][v].task_id;
		auto& task_u_id = machine[machine_id][u].task_id;
		for (int i = 0; i < _ins.get_jp()[task_v_id].size(); ++i)
		{
			auto& edge = _ins.get_jp()[task_v_id][i];
			auto& jp_v_task = _ins.get_tasks()[edge.v];
			int cc_time = D[jp_v_task.belong_machine][machine_id] * edge.w * cfg.ccr;
			//if (task_u_id == edge.v || machine[machine_id][u].R < machine[jp_v_task.belong_machine][jp_v_task.machine_seq].R + cc_time) return false;
			//if (task_u_id == edge.v || machine[machine_id][u].get_r_time() < machine[jp_v_task.belong_machine][jp_v_task.machine_seq].get_r_time() + cc_time) return false;

			if (task_u_id == edge.v || machine[machine_id][u].get_r_time() < machine[jp_v_task.belong_machine][jp_v_task.machine_seq].get_r_time()) return false;
		}
		return true;
	}

	int try_to_move_seq_front(int machine_id, int u, int v)
	{
		int length = v - u + 1;
		int u_task_id = machine[machine_id][u].task_id;
		int v_task_id = machine[machine_id][v].task_id;
		vector<int> R_copy(length), Q_copy(length);
		for (int i = u; i <= v; ++i)
		{
			R_copy[i - u] = machine[machine_id][i].R;
			Q_copy[i - u] = machine[machine_id][i].Q;
		}
		/********************* ���� R_copy *********************/
		// ���� R_copy[v]
		int jp = 0, mp = 0;
		for (int i = 0; i < _ins.get_jp()[v_task_id].size(); ++i)
		{
			auto& edge = _ins.get_jp()[v_task_id][i];
			auto& jp_v_task = _ins.get_tasks()[edge.v];
			int cc_time = D[machine_id][jp_v_task.belong_machine] * edge.w * cfg.ccr;
			jp = max(jp, machine[jp_v_task.belong_machine][jp_v_task.machine_seq].get_r_time() + cc_time);
		}
		if (u >= 1) mp = machine[machine_id][u - 1].get_r_time();
		R_copy[v - u] = max(jp, mp);

		// ���� R_copy[u]
		jp = 0, mp = 0;
		mp = R_copy[v - u] + machine[machine_id][v].dura_time;
		for (int i = 0; i < _ins.get_jp()[u_task_id].size(); ++i)
		{
			auto& edge = _ins.get_jp()[u_task_id][i];
			auto& jp_u_task = _ins.get_tasks()[edge.v];
			int cc_time = D[machine_id][jp_u_task.belong_machine] * edge.w * cfg.ccr;
			jp = max(jp, machine[jp_u_task.belong_machine][jp_u_task.machine_seq].get_r_time() + cc_time);
		}
		R_copy[u - u] = max(jp, mp);

		// w in [l1, lk]
		for (int i = u + 1; i <= v - 1; ++i)
		{
			jp = 0, mp = 0;
			mp = R_copy[i - 1 - u] + machine[machine_id][i - 1].dura_time;
			int w_task_id = machine[machine_id][i].task_id;
			for (int j = 0; j < _ins.get_jp()[w_task_id].size(); ++j)
			{
				auto& edge = _ins.get_jp()[w_task_id][j];
				auto& jp_w_task = _ins.get_tasks()[edge.v];
				int cc_time = D[machine_id][jp_w_task.belong_machine] * edge.w * cfg.ccr;
				jp = max(jp, machine[jp_w_task.belong_machine][jp_w_task.machine_seq].get_r_time() + cc_time);
			}
			R_copy[i - u] = max(jp, mp);
		}
		/*********************���� Q_copy ***********************/
		// Q_copy[lk], lk = v - 1
		int lk = v - 1;
		int lk_task_id = machine[machine_id][lk].task_id;
		int js = 0, ms = 0;
		for (int i = 0; i < _ins.get_js()[lk_task_id].size(); ++i)
		{
			auto& edge = _ins.get_js()[lk_task_id][i];
			auto& js_lk_task = _ins.get_tasks()[edge.v];
			int cc_time = D[machine_id][js_lk_task.belong_machine] * edge.w * cfg.ccr;
			js = max(js, machine[js_lk_task.belong_machine][js_lk_task.machine_seq].get_q_time() + cc_time);
		}
		if (v + 1 <= machine[machine_id].size() - 1)
		{
			ms = machine[machine_id][v + 1].get_q_time();
		}
		Q_copy[lk - u] = max(js, ms);

		// w in [u, l(k-1)]
		for (int i = lk - 1; i >= u; --i)
		{
			js = 0, ms = 0;
			ms = Q_copy[i + 1 - u] + machine[machine_id][i + 1].dura_time;
			int w_task_id = machine[machine_id][i].task_id;
			for (int j = 0; j < _ins.get_js()[w_task_id].size(); ++j)
			{
				auto& edge = _ins.get_js()[w_task_id][j];
				auto& js_w_task = _ins.get_tasks()[edge.v];
				int cc_time = D[machine_id][js_w_task.belong_machine] * edge.w * cfg.ccr;
				js = max(js, machine[js_w_task.belong_machine][js_w_task.machine_seq].get_q_time() + cc_time);
			}
			Q_copy[i - u] = max(js, ms);
		}

		// Q_copy[v] 
		js = 0, ms = 0;
		ms = machine[machine_id][u].dura_time + Q_copy[u - u];
		for (int i = 0; i < _ins.get_js()[v_task_id].size(); ++i)
		{
			auto& edge = _ins.get_js()[v_task_id][i];
			auto& js_v_task = _ins.get_tasks()[edge.v];
			int cc_time = D[js_v_task.belong_machine][machine_id] * edge.w * cfg.ccr;
			js = max(js, machine[js_v_task.belong_machine][js_v_task.machine_seq].get_q_time() + cc_time);
		}
		Q_copy[v - u] = max(js, ms);

		int new_makespan = 0;
		for (int i = u; i <= v; ++i)
		{
			new_makespan = max(new_makespan, R_copy[i - u] + Q_copy[i - u] + machine[machine_id][i].dura_time);
		}
		return new_makespan;
	}

	MoveSeq find_move_seq(int iter) // Ѱ��һ̨�����ڲ����ŵĽ���˳��
	{
		MoveSeq tabu_move, not_tabu_move;
		best_not_tabu_move_seqs.clear();
		best_tabu_move_seqs.clear();
		tabu_list.update_tabu_table(iter);
		for (int i = 0; i < critical_block.size(); ++i)
		{
			for (int j = 0; j < critical_block[i].size(); ++j)
			{
				auto& cur_cri_block = critical_block[i][j];
				int start_index = cur_cri_block.start_index;
				int end_index = cur_cri_block.end_index;

				int block_size = cur_cri_block.get_block_size();
				if (block_size == 1) continue;
				else if (block_size == 2)
				{
					if (is_move_seq_backward_legal(i, start_index, end_index))
					{
						int cur_makespan = try_to_move_seq_back(i, start_index, end_index);
						if (tabu_list.is_tabu(i, start_index, end_index, false, machine))
						{
							// ������
							if (cur_makespan < tabu_move.obj)
							{
								best_tabu_move_seqs.clear();
								tabu_move = { i, start_index, end_index, false, cur_makespan };
								best_tabu_move_seqs.emplace_back(i, start_index, end_index, false, cur_makespan);
							}
							else if (cur_makespan == tabu_move.obj)
							{
								best_tabu_move_seqs.emplace_back(i, start_index, end_index, false, cur_makespan);
							}
						}
						else
						{
							// û�б�����
							if (cur_makespan < not_tabu_move.obj)
							{
								best_not_tabu_move_seqs.clear();
								not_tabu_move = { i, start_index, end_index, false, cur_makespan };
								best_not_tabu_move_seqs.emplace_back(i, start_index, end_index, false, cur_makespan);
							}
							else if (cur_makespan == not_tabu_move.obj)
							{
								best_not_tabu_move_seqs.emplace_back(i, start_index, end_index, false, cur_makespan);
							}
						}
					}
				}
				else // �ж��Ԫ��
				{
					// ��ÿһ���м�Ԫ�� k �������ƶ�����ǰ�� start_index ������ end_index
					for (int k = start_index + 1; k <= end_index - 1; ++k)
					{
						// ���ȳ��Խ� k ����ƶ��� end_index ֮��
						if (is_move_seq_backward_legal(i, k, end_index))
						{
							int cur_makespan = try_to_move_seq_back(i, k, end_index);
							if (tabu_list.is_tabu(i, k, end_index, false, machine))
							{
								// ������
								if (cur_makespan < tabu_move.obj)
								{
									best_tabu_move_seqs.clear();
									tabu_move = { i, k, end_index, false, cur_makespan };
									best_tabu_move_seqs.emplace_back(i, k, end_index, false, cur_makespan);
								}
								else if (cur_makespan == tabu_move.obj)
								{
									best_tabu_move_seqs.emplace_back(i, k, end_index, false, cur_makespan);
								}
							}
							else
							{
								// û�б�����
								if (cur_makespan < not_tabu_move.obj)
								{
									best_not_tabu_move_seqs.clear();
									not_tabu_move = { i, k, end_index, false, cur_makespan };
									best_not_tabu_move_seqs.emplace_back(i, k, end_index, false, cur_makespan);
								}
								else if (cur_makespan == not_tabu_move.obj)
								{
									best_not_tabu_move_seqs.emplace_back(i, k, end_index, false, cur_makespan);
								}
							}
						}
						// Ȼ���Խ� k ��ǰ�ƶ��� start_index ֮ǰ
						if (is_move_seq_forward_legal(i, start_index, k))
						{
							int cur_makespan = try_to_move_seq_front(i, start_index, k);
							if (tabu_list.is_tabu(i, start_index, k, true, machine))
							{
								// ��������
								if (cur_makespan < tabu_move.obj)
								{
									best_tabu_move_seqs.clear();
									tabu_move = { i, start_index, k, true, cur_makespan };
									best_tabu_move_seqs.emplace_back(i, start_index, k, true, cur_makespan);
								}
								else if (cur_makespan == tabu_move.obj)
								{
									best_tabu_move_seqs.emplace_back(i, start_index, k, true, cur_makespan);
								}
							}
							else
							{
								// û�б�����
								if (cur_makespan < not_tabu_move.obj)
								{
									best_not_tabu_move_seqs.clear();
									not_tabu_move = { i, start_index, k, true, cur_makespan };
									best_not_tabu_move_seqs.emplace_back(i, start_index, k, true, cur_makespan);
								}
								else if (cur_makespan == not_tabu_move.obj)
								{
									best_not_tabu_move_seqs.emplace_back(i, start_index, k, true, cur_makespan);
								}
							}
						}
					}
					// ���� start_index �ƶ�����������Ԫ��֮����Ϊ֮ǰ���Թ��� strat_index + 1 �ƶ�����ǰ����
						// ���Կ���ֱ�ӽ� start_index �ƶ��� start_index + 2 �����λ����
					for (int k = start_index + 2; k <= end_index; ++k)
					{
						if (is_move_seq_backward_legal(i, start_index, k))
						{
							int cur_makespan = try_to_move_seq_back(i, start_index, k);
							if (tabu_list.is_tabu(i, start_index, k, false, machine))
							{
								// ������
								if (cur_makespan < tabu_move.obj)
								{
									best_tabu_move_seqs.clear();
									tabu_move = { i, start_index, k, false, cur_makespan };
									best_tabu_move_seqs.emplace_back(i, start_index, k, false, cur_makespan);
								}
								else if (cur_makespan == tabu_move.obj)
								{
									best_tabu_move_seqs.emplace_back(i, start_index, k, false, cur_makespan);
								}
							}
							else
							{
								// û�б�����
								if (cur_makespan < not_tabu_move.obj)
								{
									best_not_tabu_move_seqs.clear();
									not_tabu_move = { i, start_index, k, false, cur_makespan };
									best_not_tabu_move_seqs.emplace_back(i, start_index, k, false, cur_makespan);
								}
								else if (cur_makespan == not_tabu_move.obj)
								{
									best_not_tabu_move_seqs.emplace_back(i, start_index, k, false, cur_makespan);
								}
							}
						}
						else
						{
							break;
						}
					}
					// ���� end_index �ƶ�����������Ԫ��֮ǰ����Ϊ֮ǰ���Թ��� end_index - 1 �ƶ���������
						// ���Կ���ֱ�ӽ� end_index �ƶ��� end_index - 2 ǰ���λ����
					for (int k = end_index - 2; k >= start_index; --k)
					{
						if (is_move_seq_forward_legal(i, k, end_index))
						{
							int cur_makespan = try_to_move_seq_front(i, k, end_index);
							if (tabu_list.is_tabu(i, k, end_index, true, machine))
							{
								// ��������
								if (cur_makespan < tabu_move.obj)
								{
									best_tabu_move_seqs.clear();
									tabu_move = { i, k, end_index, true, cur_makespan };
									best_tabu_move_seqs.emplace_back(i, k, end_index, true, cur_makespan);
								}
								else if (cur_makespan == tabu_move.obj)
								{
									best_tabu_move_seqs.emplace_back(i, k, end_index, true, cur_makespan);
								}
							}
							else
							{
								// û�б�����
								if (cur_makespan < not_tabu_move.obj)
								{
									best_not_tabu_move_seqs.clear();
									not_tabu_move = { i, k, end_index, true, cur_makespan };
									best_not_tabu_move_seqs.emplace_back(i, k, end_index, true, cur_makespan);
								}
								else if (cur_makespan == not_tabu_move.obj)
								{
									best_not_tabu_move_seqs.emplace_back(i, k, end_index, true, cur_makespan);
								}
							}
						}
						else
						{
							break;
						}
					}
				}
			}
		}

		if (best_tabu_move_seqs.empty() && best_not_tabu_move_seqs.empty())
		{
			tabu_move = { -1, -1, -1, false, -1 };
			return tabu_move;
		}
		if ((tabu_move.obj < not_tabu_move.obj && tabu_move.obj < best_makespan) || best_not_tabu_move_seqs.empty())
		{
			// �������
			auto& select_move = best_tabu_move_seqs[fast_rand(best_tabu_move_seqs.size())];
			tabu_list.remove_tabu(select_move.machine_id, select_move.u, select_move.v, select_move.is_forward, machine);
			tabu_list.make_tabu(select_move.machine_id, select_move.u, select_move.v, iter, fast_rand(tabu_list.get_L(), tabu_list.get_R()), machine);
			return select_move;
		}
		else
		{
			// ѡ��û�б����ɵ�
			auto& select_move = best_not_tabu_move_seqs[fast_rand(best_not_tabu_move_seqs.size())];
			tabu_list.make_tabu(select_move.machine_id, select_move.u, select_move.v, iter, fast_rand(tabu_list.get_L(), tabu_list.get_R()), machine);
			return select_move;
		}
	}

	// �� v ��ǰ�����
	int score_move_machine_to_loc(int m1, int u, int m2, int v, int jp, int js)
	{
		int mp, ms;
		if (v == 0) mp = 0;
		else mp = machine[m2][v - 1].get_r_time();

		if (v == machine[m2].size()) ms = 0;
		else ms = machine[m2][v].get_q_time();

		int R_copy = max(jp, mp);
		int Q_copy = max(js, ms);
		return max(machine[m1][u].dura_time + R_copy + Q_copy, machine[m1][u].R + machine[m1][u].Q);
	}

	// TODO:: �����費��Ҫ�� Start �� End �ų���
	// ��Ϊ���� -1 �ˣ�����ֵ��Ҫ���� v �ĺ���
	MoveMachine try_to_move_machine(int m1, int u, int m2)
	{
		int L = 0, R = machine[m2].size();
		if (machine[m2].size() >= 1)
		{
			if (machine[m2][0].task_id == Start) L++;
			if (machine[m2].back().task_id == End) R--;
		}
		int left = L, right = R - 1;
		int jp = 0, js = 0;
		MoveMachine ret;
		// ��һ�� JP(u)
		// �����˵� ���� jp(u).Q - cc_time >= left.Q, ������left��ǰ�����
		// ��֤�� JP(u) ���֮�����ݴ��͵� u 
		//int min_jp_sub_cc_time = INF;
		int min_jp_sub_cc_time = 0;
		int u_task_id = machine[m1][u].task_id;
		set<int> jp_task_ids;
		for (int i = 0; i < _ins.get_jp()[u_task_id].size(); ++i)
		{
			auto& edge = _ins.get_jp()[u_task_id][i];
			auto& jp_u_task = _ins.get_tasks()[edge.v];
			jp_task_ids.insert(edge.v);
			int cc_time = D[m2][jp_u_task.belong_machine] * edge.w * cfg.ccr;
			//min_jp_sub_cc_time = min(min_jp_sub_cc_time, machine[jp_u_task.belong_machine][jp_u_task.machine_seq].Q - cc_time);
			min_jp_sub_cc_time = max(min_jp_sub_cc_time, machine[jp_u_task.belong_machine][jp_u_task.machine_seq].get_r_time());
			jp = max(jp, machine[jp_u_task.belong_machine][jp_u_task.machine_seq].get_r_time() + cc_time);
		}
		// rigth u left
		// ֮��Ķ����� jp(u).Q - cc_time >= left.Q
		// ��� left �� Start �Ļ���Start.Q �����ģ��϶������㣬++֮��Ͳ��õ�����
		//while (left < machine[m2].size() && ((machine[m2][left].Q > min_jp_sub_cc_time || jp_task_ids.count(machine[m2][left].task_id))))
		//{
		//	left++;
		//}
		
		//while (left < R && ((machine[m2][left].Q > min_jp_sub_cc_time || jp_task_ids.count(machine[m2][left].task_id))))
		//{
		//	left++;
		//}

		while (left < R && ((machine[m2][left].get_r_time() < min_jp_sub_cc_time || jp_task_ids.count(machine[m2][left].task_id))))
		{
			left++;
		}

		// ��һ�� JS(u)
		// ����Ҷ˵㣬���� right.Q >= js(u).Q + cc_time; �����Ļ������� right �ĺ�����룬
		// �� right + 1 ��ǰ�����
		// ��֤�� u ���֮���ܹ������ݴ��͵� js(u)
		jp_task_ids.clear();
		int max_js_add_cc_time = 0;
		auto& js_task_ids = jp_task_ids;
		for (int i = 0; i < _ins.get_js()[u_task_id].size(); ++i)
		{
			auto& edge = _ins.get_js()[u_task_id][i];
			auto& js_u_task = _ins.get_tasks()[edge.v];
			js_task_ids.insert(edge.v);
			int cc_time = D[m2][js_u_task.belong_machine] * edge.w * cfg.ccr;
			//max_js_add_cc_time = max(max_js_add_cc_time, machine[js_u_task.belong_machine][js_u_task.machine_seq].Q + cc_time);
			max_js_add_cc_time = max(max_js_add_cc_time, machine[js_u_task.belong_machine][js_u_task.machine_seq].get_q_time());
			js = max(js, machine[js_u_task.belong_machine][js_u_task.machine_seq].get_q_time() + cc_time);
		}

		// right u left
		// ֮ǰ��ȫ������ right.Q >= js(u).Q
		// 
		//while (right >= 0 && ((machine[m2][right].Q < max_js_add_cc_time) || js_task_ids.count(machine[m2][right].task_id)))
		//{
		//	right--;
		//}
		
		//while (right >= L && ((machine[m2][right].Q < max_js_add_cc_time) || js_task_ids.count(machine[m2][right].task_id)))
		//{
		//	right--;
		//}

		while (right >= L && ((machine[m2][right].get_q_time() < max_js_add_cc_time) || js_task_ids.count(machine[m2][right].task_id)))
		{
			right--;
		}
		right++; // ��ΪҪ��right֮����룬Ϊ������ͳһ�ԣ�������һ�£������Ļ���ȫ��ǰ����
		// [left, right]
		if (right - left + 1 <= 0)
		{
			// �����ںϷ��ĸı�
			ret.obj = -1;
			return ret;
		}

		best_move_machine_score_eqs.clear();
		for (int i = left; i <= right; ++i)
		{
			int obj = score_move_machine_to_loc(m1, u, m2, i, jp, js);
			if (obj < ret.obj)
			{
				ret = { m1, u, m2, i - 1, obj };
				best_move_machine_score_eqs.clear();
				best_move_machine_score_eqs.emplace_back(m1, u, m2, i - 1, obj);
			}
			else if (obj == ret.obj)
			{
				best_move_machine_score_eqs.emplace_back(m1, u, m2, i - 1, obj);
			}
		}
		return best_move_machine_score_eqs[fast_rand(best_move_machine_score_eqs.size())];
	}

	MoveMachine find_move_machine(int iter)
	{
		MoveMachine tabu_move, not_tabu_move;
		best_tabu_move_machines.clear();
		best_not_tabu_move_machines.clear();
		for (int i = 0; i < critical_block.size(); ++i)
		{
			for (int j = 0; j < critical_block[i].size(); ++j)
			{
				int start_index = critical_block[i][j].start_index;
				int end_index = critical_block[i][j].end_index;
				for (int k = start_index; k <= end_index; ++k)
				{
					auto& k_task_id = machine[i][k].task_id;
					auto& k_task = _ins.get_tasks()[k_task_id];
					for (int l = 0; l < machine.size(); ++l)
					{
						// ���Խ� k_task ת�Ƶ� L �Ż�����
						if (l == k_task.belong_machine) continue;
						MoveMachine temp = try_to_move_machine(i, k, l);
						if (temp.v == -1 && temp.obj == -1) continue;
						if (!tabu_machine.is_tabu(k_task_id, l, iter))
						{
							// û�б�����
							if (temp.obj < not_tabu_move.obj)
							{
								not_tabu_move = temp;
								best_not_tabu_move_machines.clear();
								best_not_tabu_move_machines.push_back(temp);
							}
							else if (temp.obj == not_tabu_move.obj)
							{
								best_not_tabu_move_machines.push_back(temp);
							}
						}
						else
						{
							// ������
							if (temp.obj < tabu_move.obj)
							{
								tabu_move = temp;
								best_tabu_move_machines.clear();
								best_tabu_move_machines.push_back(temp);
							}
							else if (temp.obj == tabu_move.obj)
							{
								best_tabu_move_machines.push_back(temp);
							}
						}
					}
				}
			}
		}

		if (best_not_tabu_move_machines.empty() && best_tabu_move_machines.empty()) { not_tabu_move.obj = -1; return not_tabu_move; }
		if (tabu_move.obj < not_tabu_move.obj && tabu_move.obj < best_makespan || best_not_tabu_move_machines.empty())
		{
			// �����Ϊ������һ��������ֱ��ʹ��tabu_move����
			tabu_move = best_tabu_move_machines[fast_rand(best_tabu_move_machines.size())];
			int u_task_id = machine[tabu_move.m1][tabu_move.u].task_id;
			tabu_machine.make_tabu(u_task_id, tabu_move.m1, iter, fast_rand(tabu_machine.get_L(), tabu_machine.get_R()));
			return tabu_move;
		}
		else
		{
			// �������ɵ�
			not_tabu_move = best_not_tabu_move_machines[fast_rand(best_not_tabu_move_machines.size())];
			int u_task_id = machine[not_tabu_move.m1][not_tabu_move.u].task_id;
			tabu_machine.make_tabu(u_task_id, not_tabu_move.m1, iter, fast_rand(tabu_machine.get_L(), tabu_machine.get_R()));
			return not_tabu_move;
		}
	}

	void make_move_seq(MoveSeq move_seq)
	{
		int task_u_id = machine[move_seq.machine_id][move_seq.u].task_id;
		int task_v_id = machine[move_seq.machine_id][move_seq.v].task_id;
		if (move_seq.is_forward)
		{	// �� v Ų�� u ��ǰ�棬����Ҫ�޸� task ����
			_ins.get_tasks()[task_v_id].machine_seq = move_seq.u;
			machine[move_seq.machine_id].insert(machine[move_seq.machine_id].begin() + move_seq.u, machine[move_seq.machine_id][move_seq.v]);
			machine[move_seq.machine_id].erase(machine[move_seq.machine_id].begin() + move_seq.v + 1);
			// ɾ����Ӱ��Ԫ�ص� machine_seq;
			// [u + 1, v] seq++
			for (int i = move_seq.u + 1; i <= move_seq.v; ++i)
			{
				int task_id = machine[move_seq.machine_id][i].task_id;
				_ins.get_tasks()[task_id].machine_seq = i;
			}
		}
		else
		{
			// �� u Ų�� v �ĺ��棬����Ҫ�޸� task ����
			_ins.get_tasks()[task_u_id].machine_seq = move_seq.v;
			machine[move_seq.machine_id].insert(machine[move_seq.machine_id].begin() + move_seq.v + 1, machine[move_seq.machine_id][move_seq.u]);
			machine[move_seq.machine_id].erase(machine[move_seq.machine_id].begin() + move_seq.u);
			// [u, v - 1] seq--
			for (int i = move_seq.u; i <= move_seq.v - 1; ++i)
			{
				int task_id = machine[move_seq.machine_id][i].task_id;
				_ins.get_tasks()[task_id].machine_seq = i;
			}
		}
	}

	void make_move_machine(MoveMachine move_machine)
	{
		//�ӻ���m1��uλ��Ų��������m2��v + 1λ��,ͬʱ�޸ĸ���ҵѡ��Ļ���ֵ
		int u_task_id = machine[move_machine.m1][move_machine.u].task_id;
		_ins.get_tasks()[u_task_id].belong_machine = move_machine.m2;
		_ins.get_tasks()[u_task_id].machine_seq = move_machine.v + 1;
		machine[move_machine.m2].insert(machine[move_machine.m2].begin() + (move_machine.v + 1), machine[move_machine.m1][move_machine.u]);
		// m2, [v + 2, end] seq++
		for (int i = move_machine.v + 2; i < machine[move_machine.m2].size(); ++i)
		{
			int task_id = machine[move_machine.m2][i].task_id;
			_ins.get_tasks()[task_id].machine_seq = i;
		}
		// m1 [u, end] seq--
		machine[move_machine.m1].erase(machine[move_machine.m1].begin() + move_machine.u);
		for (int i = move_machine.u; i < machine[move_machine.m1].size(); ++i)
		{
			int task_id = machine[move_machine.m1][i].task_id;
			_ins.get_tasks()[task_id].machine_seq = i;
		}
	}

	bool tabu_search()
	{
		int already_fail_times = 0;
		int iter = 1;
		tabu_list.init();
		tabu_machine.init();
		bool is_update_best_makespan = false;
		while ((already_fail_times < cfg.max_tabu_times) && !timer.isTimeout())
		{
			MoveSeq move_seq = find_move_seq(iter);
			MoveMachine move_machine = find_move_machine(iter);
			
			
			if (move_seq.obj != -1 && move_machine.obj != -1)
			{
				if (move_seq.obj < move_machine.obj)
				{
					watch_move_seq(move_seq);
					make_move_seq(move_seq);
				}
				else if (move_seq.obj > move_machine.obj)
				{
					watch_move_machine(move_machine);
					make_move_machine(move_machine);
				}
				else if (move_seq.obj == move_machine.obj)
				{
					if (rand() % 2 == 0)
					{
						watch_move_seq(move_seq);
						make_move_seq(move_seq);
					}
					else
					{
						watch_move_machine(move_machine);
						make_move_machine(move_machine);
					}
				}
			}
			else if (move_seq.obj == -1 && move_machine.obj == -1)
			{
				watch("�Ҳ�����");
				break;
			}
			else if (move_seq.obj == -1 && move_machine.obj != -1)
			{
				watch_move_machine(move_machine);
				make_move_machine(move_machine);
			}
			else if (move_seq.obj != -1 && move_machine.obj == -1)
			{
				watch_move_seq(move_seq);
				make_move_seq(move_seq);
			}

			sum_critical_path();
			if (check_obj(_all_iteration))
			{
				is_update_best_makespan = true;
				already_fail_times = 0;
			}
			else
			{
				already_fail_times++;
			}
			iter++;
			_all_iteration++;
		}
		return is_update_best_makespan;
	}

	void test_R()
	{
		
		for (int i = 0; i < machine.size(); ++i)
		{
			int p = -1;
			for (int j = 0; j < machine[i].size(); ++j)
			{
				assert(machine[i][j].R >= p);
				p = machine[i][j].R;
			}
		}
	}

	// [left = -1, right = -1]   size = 0
	// [left = randloc, right = left + 1]�� left in [0, size - 2]; right in [1, size - 1]
	// [left = size - 1, right = -1]
	// = -1 ��ʱ�����ж��ˣ�ֱ�ӺϷ�
	bool is_move_machine_legal(int m1, int u, int m2, int left, int right)
	{
		int u_task_id = machine[m1][u].task_id;
		if (left == -1)
		{
			// ֻ�� >= 1 ��ʱ��Start�ſ��ܻ�����ڿ���
			if (machine[m2].size() >= 1 && machine[m2][0].task_id == Start) return false;
		}
		else if (left <= machine[m2].size() - 1)
		{
			if (left == machine[m2].size() - 1)
			{
				if (machine[m2].back().task_id == End) return false;
			}
			// js_u
			//int max_js_q = 0;
			//set<int> js_u_ids;
			//for (int i = 0; i < _ins.get_js()[u_task_id].size(); ++i)
			//{
			//	auto& edge = _ins.get_js()[u_task_id][i];
			//	js_u_ids.insert(edge.v);
			//	auto& js_u_task = _ins.get_tasks()[edge.v];
			//	int cc_time = D[m1][js_u_task.belong_machine] * edge.w * cfg.ccr;
			//	max_js_q = max(max_js_q, machine[js_u_task.belong_machine][js_u_task.machine_seq].Q + cc_time);
			//}
			//if ((max_js_q > machine[m2][left].Q) || js_u_ids.count(machine[m2][left].task_id))
			//	return false;

			for (int i = 0; i < _ins.get_js()[u_task_id].size(); ++i)
			{
				auto& edge = _ins.get_js()[u_task_id][i];
				auto& js_u_task = _ins.get_tasks()[edge.v];
				int cc_time = D[m1][js_u_task.belong_machine] * edge.w * cfg.ccr;
				//if (machine[js_u_task.belong_machine][js_u_task.machine_seq].Q + cc_time > machine[m2][left].Q
				//	|| edge.v == machine[m2][left].task_id) return false;

				if (machine[js_u_task.belong_machine][js_u_task.machine_seq].get_q_time() > machine[m2][left].get_q_time()
					|| edge.v == machine[m2][left].task_id) return false;
			}
			
		}
		//else if (left == machine[m2].size() - 1)
		//{
		//	if (machine[m2].back().task_id == End) return false;
		//}

		if (right != -1)
		{
			assert(machine[m1][u].task_id != End);
			assert(machine[m1][u].task_id != Start);

			/*int min_jp_q_cc_time = INF;
			set<int> jp_u_ids;
			for (int i = 0; i < _ins.get_jp()[u_task_id].size(); ++i)
			{
				auto& edge = _ins.get_jp()[u_task_id][i];
				jp_u_ids.insert(edge.v);
				auto& jp_u_task = _ins.get_tasks()[edge.v];
				int cc_time = D[m1][jp_u_task.belong_machine] * edge.w * cfg.ccr;
				min_jp_q_cc_time = min(min_jp_q_cc_time, machine[jp_u_task.belong_machine][jp_u_task.machine_seq].Q - cc_time);
			}
			if ((min_jp_q_cc_time < machine[m2][right].Q) || jp_u_ids.count(machine[m1][u].task_id))
				return false;*/

			for (int i = 0; i < _ins.get_jp()[u_task_id].size(); ++i)
			{
				auto& edge = _ins.get_jp()[u_task_id][i];
				auto& jp_u_task = _ins.get_tasks()[edge.v];
				int cc_time = D[m1][jp_u_task.belong_machine] * edge.w * cfg.ccr;
				//if (machine[jp_u_task.belong_machine][jp_u_task.machine_seq].Q - cc_time < machine[m2][right].Q
				//	|| edge.v == machine[m2][right].task_id) return false;

				if (machine[jp_u_task.belong_machine][jp_u_task.machine_seq].get_r_time() > machine[m2][right].get_r_time()
					|| edge.v == machine[m2][right].task_id) return false;

			}	
		}
		return true;
	}

	void randmove()
	{
		MoveSeq move_seq;
		MoveMachine move_machine;
		// ���һ������������ Move_seq
		bool flag = true;
		for (int i = 0; i < machine.size(); ++i)
		{
			int num = machine[i].size();
			if (num <= 1) continue;
			if (machine[i][0].task_id == Start) num--;
			if (machine[i].back().task_id == End) num--;
			if (num >= 2)
			{
				flag = false;
				break;
			}
		}
		if (flag) return;
		// ���Խ��� move_seq������ 50 ��
		bool find_move_seq = false, find_move_machine = false;
		int x = 0;
		while (x++ < cfg.randmove_times)
		{
			int select_machine = -1;
			while (true)
			{
				select_machine = fast_rand(machine.size());
				if (machine[select_machine].size() <= 1) continue;
				int num = machine[select_machine].size();
				if (machine[select_machine][0].task_id == Start) num--;
				if (machine[select_machine].back().task_id == End) num--;
				if (num >= 2) break;
			}
			int l = 0, r = machine[select_machine].size();
			if (machine[select_machine][0].task_id == Start) l++;
			if (machine[select_machine].back().task_id == End) r--;
			int u = fast_rand(l, r);
			int v = fast_rand(l, r);
			while(u == v) v = fast_rand(l, r);
			if (u > v) swap(u, v);
			bool back = is_move_seq_backward_legal(select_machine, u, v);
			bool front = is_move_seq_forward_legal(select_machine, u, v);
			if (back || front)
			{
				int num = 1;
				if (back)
				{
					num++;
					move_seq = MoveSeq(select_machine, u, v, false, 0);
				}
				if (front)
				{
					if (fast_rand(num) == 0)
						move_seq = MoveSeq(select_machine, u, v, true, 0);
				}
				find_move_seq = true;
				break;	// ֻҪ����һ���Ϸ��ͻ�ֱ�� break ��
			}
		}
		if(move_seq.u == -1 && move_seq.v == -1) watch("����ƶ��쳣��");
		if (fast_rand(2) == 0 && find_move_seq)
		{
			//cerr << "ִ���� make_move_seq \t";
			//watch(machine[move_seq.machine_id].size(), move_seq.machine_id, move_seq.u, move_seq.v, move_seq.is_forward);
			make_move_seq(move_seq);
		}
		else
		{
			int x = 0;
			while (x < cfg.randmove_times) // ���� 50 ��move_machine��ʧ�ܵĻ��� move_seq
			{
				int m1 = -1;
				while (true)
				{
					m1 = fast_rand(machine.size());
					if (machine[m1].size() == 0) continue;
					int num = machine[m1].size();
					if (machine[m1][0].task_id == Start) num--;
					if (machine[m1].back().task_id == End) num--;
					if (num >= 1) break;
				}
				int l = 0, r = machine[m1].size();
				if (machine[m1][0].task_id == Start) l++;
				if (machine[m1].back().task_id == End) r--;
				int u = fast_rand(l, r);	// ��֤�� u ������ Start �� End
				int m2 = fast_rand(machine.size());
				while (m2 == m1) m2 = fast_rand(machine.size());

				int left, right;
				if (machine[m2].size() == 0)
				{
					left = right = -1;		// return true;
				}
				else
				{
					// [0, size]����Ϊ�����λ�û��Ԫ�ظ�����һ��
					int randloc = fast_rand(machine[m2].size() + 1);
					// [-1, size - 1]	��Ϊ�Ǻ���룬���Բ��õ���
					randloc--;
					left = randloc;
					if (left == machine[m2].size() - 1) right = -1;
					else right = left + 1;
				}
				// [left = -1, right = -1]
				// [left = randloc, right = left + 1]�� left in [-1, size - 1]
				// [left = size - 1, right = -1]
				// = -1 ��ʱ�����ж��ˣ�ֱ�ӺϷ�
				// Ҫ��֤�����λ�ò��� Start ֮ǰ��End ֮��
				if (is_move_machine_legal(m1, u, m2, left, right))
				{
					move_machine = MoveMachine(m1, u, m2, left, 0);
					//cerr << "ִ���� make_move_machine \t";
					//watch(machine[move_machine.m1].size(), machine[move_machine.m2].size(), move_machine.m1, move_machine.u, move_machine.m2, move_machine.v);
					make_move_machine(move_machine);
					break;
				}
				x++;
			}
			if (move_machine.m1 == -1 && move_machine.v == -1) watch("����ƶ�---����---�쳣��");
			if (x == cfg.randmove_times && find_move_seq) make_move_seq(move_seq);
		}
	}

	void run() {
		init_solution();
		sum_critical_path();
		check_obj(_all_iteration);
		int pre_makespan = makespan;
		
		while (!timer.isTimeout() && best_makespan != _ins.get_best_ans()) 
		{
			tabu_search();
			watch(timer.getDuration(), pre_makespan, makespan, best_makespan);
			perturbation(cfg.randmove_max_iter);
			pre_makespan = makespan;
		}
	}

	void perturbation(int times)
	{
		//machine = best_machine;
		//for (int i = 0; i < _ins.get_vir_task_num(); ++i)
		//{
		//	_ins.get_tasks()[i].belong_machine = task_of_machine[i].first;
		//	_ins.get_tasks()[i].machine_seq = task_of_machine[i].second;
		//}
		//sum_critical_path();
		for (int i = 0; i < times; ++i)
		{
			randmove();
			sum_critical_path();
		}
	}

	// �����д���ļ�
	void record_sol(const string& sol_path) const {
		freopen(sol_path.c_str(), "w", stdout);
		record_sol();
	}
	// �����д�����̨
	void record_sol() const
	{
		cout << best_makespan << endl;
		for (int i = 0; i < best_machine.size(); ++i)
		{
			for (const auto& item : best_machine[i])
			{
				if (item.task_id == Start || item.task_id == End) continue;
				cout << item.task_id << " " << item.R << " " << item.dura_time << "    ";
			}
			cout << endl;
		}
	}
	// ���ӻ����
	void draw_sol(const string& html_path) const {

	}
	//// ���ӻ����
	//void draw_sol(const string& html_path) const {
	//	utils_visualize_drawer::Drawer html_drawer(html_path, _width, _height);
	//	for (auto& dst_node : _dst) {
	//		html_drawer.rect(dst_node.lb_point.x, dst_node.lb_point.y,
	//			dst_node.width, dst_node.height, dst_node.rotation == Rotation::_90_);
	//	}
	//}

	// д���׼������
	void record_log() {
		// д���ͷ
		writeLog(LineType::HEAD, "date", "instance", "makespan", "cfg.ccr", "total_iterTimes", "optimal_iterTimes", "randomSeed", "total_time", "optimal_time");
		// д��������
		writeLog(LineType::LINE, utils::Date::to_format_str(), _env.instance_name(), best_makespan, cfg.ccr, _all_iteration, _iteration, _cfg.random_seed, timer.getDuration(), _duration);
	}
	// д����־�ļ�
	void record_log(string log_path)
	{
		// д���ͷ
		writeLog(log_path, LineType::HEAD, "date", "instance", "makespan", "cfg.ccr", "total_iterTimes", "optimal_iterTimes", "randomSeed", "total_time", "optimal_time");
		// д��������
		writeLog(log_path, LineType::LINE, utils::Date::to_format_str(), _env.instance_name(), best_makespan, cfg.ccr, _all_iteration, _iteration, _cfg.random_seed, timer.getDuration(), _duration);
	}

public:

	/// ������Ƿ���ȷ
	void check() {
		// best_machine
		queue<int> q;
		vector<int> dis(_ins.get_vir_task_num());
		vector<int> indeg(_ins.get_vir_task_num());
		for (int i = 0; i < _ins.get_vir_task_num(); ++i) indeg[i] = _ins.get_jp()[i].size();
		for (int i = 0; i < best_machine.size(); ++i)
		{
			for (int j = 0; j < best_machine[i].size(); ++j)
			{
				int task_id = best_machine[i][j].task_id;
				indeg[task_id] += (j != 0);
			}
		}
		// ���� ��һ�����е�Ŀ�ʼʱ��
		dis[Start] = 0;
		q.push(Start);
		while (!q.empty())
		{
			int u = q.front(); q.pop();
			// machine ����
			auto& task_u = _ins.get_tasks()[u];
			if (task_of_machine[u].second + 1 < best_machine[task_of_machine[u].first].size())
			{
				int ms_id = best_machine[task_of_machine[u].first][task_of_machine[u].second + 1].task_id;
				dis[ms_id] = max(dis[ms_id], dis[u] + task_u.length);
				if (--indeg[ms_id] == 0)
				{
					q.push(ms_id);
				}
			}
			for (int i = 0; i < _ins.get_js()[u].size(); ++i)
			{
				auto& edge = _ins.get_js()[u][i];
				auto& js_u_task = _ins.get_tasks()[edge.v];
				int cc_time = D[task_of_machine[u].first][task_of_machine[edge.v].first] * edge.w * cfg.ccr;
				dis[edge.v] = max(dis[edge.v], dis[u] + task_u.length + cc_time);
				if (--indeg[edge.v] == 0)
				{
					q.push(edge.v);
				}
			}
		}

		assert(dis[End] == best_makespan);
		for (int i = 0; i < _ins.get_vir_task_num(); ++i)
		{
			//auto& task_u = _ins.get_tasks()[i];
			for (int j = 0; j < _ins.get_jp()[i].size(); ++j)
			{
				auto& edge = _ins.get_jp()[i][j];
				auto& jp_u_task = _ins.get_tasks()[edge.v];
				int cc_time = D[task_of_machine[i].first][task_of_machine[edge.v].first] * edge.w * cfg.ccr;
				// check ����ǰ�����֮�󣬸�������ܿ�ʼ
				assert(dis[edge.v] + cc_time + jp_u_task.length <= dis[edge.u]);
			}
		}
		for (int i = 0; i < best_machine.size(); ++i)
		{
			int pre = 0;
			for (int j = 0; j < best_machine[i].size(); ++j)
			{
				assert(dis[best_machine[i][j].task_id] >= pre);
				pre = dis[best_machine[i][j].task_id] + best_machine[i][j].dura_time;
			}
		}
		
	}

	void watch_move_seq(MoveSeq move_seq)
	{
#ifdef SHOW_MOVE
		int u = machine[move_seq.machine_id][move_seq.u].task_id;
		int v = machine[move_seq.machine_id][move_seq.v].task_id;
		bool is_forward = move_seq.is_forward;
		watch(move_seq.obj, u, v, is_forward, move_seq.u, move_seq.v);
#endif // SHOW_MOVE
	}

	void watch_move_machine(MoveMachine move_machine)
	{
#ifdef SHOW_MOVE
		int m1 = move_machine.m1;
		int u = machine[move_machine.m1][move_machine.u].task_id;
		int m2 = move_machine.m2;
		
		int v = machine[move_machine.m2][max(move_machine.v, 0)].task_id;
		watch(move_machine.obj, m1, u, m2, v);
#endif // SHOW_MOVE
	}

};
#endif // !_SRC_SOLVE_H_