#ifndef _SRC_STRUCT_H_
#define _SRC_STRUCT_H_
#include <cassert>
/*
* 该文件，编写一些需要用到的数据结构，数据类型等
*/
struct Edge
{
	int u, v, w;
	Edge(int uu, int vv, int ww): u(uu), v(vv), w(ww){}
};

struct Task
{
	int length;
	int belong_machine;
	int machine_seq;
	int start_time;
	Task() 
	{
		length = 0;
		belong_machine = -1;
		machine_seq = -1;
		start_time = 0;
	}
	Task(int len) : length(len), belong_machine(-1), start_time(0), machine_seq(0) {}
	int get_total_time()
	{
		return this->start_time + this->length;
	}
};

struct Operation
{
	int task_id;//工件
	int dura_time;//持续时间
	int Q;
	int R;//从起点到当前操作的时间
	Operation() {}
	Operation(int id, int time, int r) :task_id(id), dura_time(time), R(r)
	{
	}
	int get_r_time()
	{
		return R + dura_time;
	}
	int get_q_time()
	{
		return Q + dura_time;
	}
	int get_r_q_time()
	{
		return R + Q + dura_time;
	}
};

struct Block
{
	int start_index, end_index;
	int get_block_size() { return end_index - start_index + 1; }
};

// N7邻域，move操作
struct MoveSeq
{
	int machine_id;
	int u;
	int v;
	bool is_forward; // 是否是向前插入
	int obj;	// 评分，分越小越好
	MoveSeq(int id, int uu, int vv, bool forward, int obj): machine_id(id), u(uu),
		v(vv), is_forward(forward), obj(obj){}
	MoveSeq()
	{
		machine_id = -1;
		u = -1;
		v = -1;
		is_forward = true;
		obj = INF;
	}
};

// 将机器m1上的u挪到m2机器v的后面
struct MoveMachine
{
	int m1, u, m2, v, obj;
	MoveMachine() 
	{
		m1 = u = m2 = v = -1;
		obj = INF;
	}
	MoveMachine(int m1, int u, int m2, int v, int obj) :
		m1(m1), u(u), m2(m2), v(v), obj(obj) {}
};

struct TabuBlock
{
	int tabu_time;
	vector<int> tasks;
	TabuBlock() { tabu_time = 0; }
};

struct TabuTable
{
private:
	int length_min, length_max;
	vector<vector<TabuBlock>> tabu_list;
public:
	
	TabuTable(){}
	TabuTable(int machine_num, int l_min, int l_max): tabu_list(machine_num), length_min(l_min), length_max(l_max) {}
	void init()
	{
		for (int i = 0; i < tabu_list.size(); ++i)
		{
			tabu_list[i].clear();
		}
	}
	void make_tabu(int machine_id, int u, int v, int iter, int len, vector<vector<Operation>>& machine)
	{
		// 禁忌当前状态
		TabuBlock tabu_block;
		for (int i = 0; i < v - u + 1; ++i)
		{
			tabu_block.tasks.push_back(machine[machine_id][u + i].task_id);
		}
		tabu_block.tabu_time = iter + len; // len 禁忌步长，需要在外面使用随机数，row出来，然后传进来。
		tabu_list[machine_id].push_back(tabu_block);
	}
	bool is_tabu(int machine_id, int u, int v, bool is_forward, vector<vector<Operation>>& machine)
	{
		return is_exist(machine_id, u, v, is_forward, machine) != -1;
	}

	// 返回该段禁忌序列的下标
	int is_exist(int machine_id, int u, int v, bool is_forward, vector<vector<Operation>>& machine)
	{
		int len = v - u + 1;
		vector<int> w(len);
		if (is_forward)
		{
			w[0] = v;
			for (int i = 0; i <= len - 2; ++i) { w[i + 1] = u + i; }
		}
		else
		{
			w[v - u] = u;
			for (int i = 0; i <= len - 2; ++i) { w[i] = u + i + 1; }
		}
		for (int i = 0; i < tabu_list[machine_id].size(); ++i)
		{
			auto& tabu_block = tabu_list[machine_id][i];
			if (tabu_block.tasks.size() >= len)
			{
				// 两个序列进行匹配，其实可以使用 kmp
				for (int s = 0; s <= tabu_block.tasks.size() - len; ++s)
				{
					bool find_flag = true;
					for (int k = 0; k < len; ++k)
					{
						if (tabu_block.tasks[s + k] != machine[machine_id][w[k]].task_id)
						{
							find_flag = false;
							break;
						}
					}
					if (find_flag)
					{
						return i;
					}
				}
			}
		}
		return -1;
	}

	void update_tabu_table(int time)
	{
		for (int i = 0; i < tabu_list.size(); ++i)
		{
			for (auto iter = tabu_list[i].begin(); iter != tabu_list[i].end();)
			{
				if (iter->tabu_time < time)
				{
					iter = tabu_list[i].erase(iter);
				}
				else { ++iter; }
			}
		}
	}
	void remove_tabu(int machine_id, int u, int v, bool is_forward, vector<vector<Operation>>& machine)
	{
		int exist_index = is_exist(machine_id, u, v, is_forward, machine);
		assert(exist_index != -1);
		tabu_list[machine_id].erase(tabu_list[machine_id].begin() + exist_index);
	}

	int get_L() { return length_min; }
	int get_R() { return length_max; }
};

struct TabuArr
{
private:
	int length_min, length_max;
	int _rows, _cols;
	vector<vector<int>> tabu_list;
public:
	int get_L() { return length_min; }
	int get_R() { return length_max; }
	TabuArr() {}
	TabuArr(int rows, int cols, int len_min, int len_max) : length_min(len_min), length_max(len_max), _rows(rows + 5), _cols(cols + 5)
	{
		tabu_list.resize(_rows);
		for (int i = 0; i < _rows; ++i) { tabu_list[i].resize(_cols); }
	}
	bool is_tabu(int x, int y, int iter) { return tabu_list[x][y] >= iter; }
	void make_tabu(int x, int y, int iter, int len) { tabu_list[x][y] = iter + len; }
	void remove_tabu(int x, int y) { tabu_list[x][y] = 0; }
	void init() { for (int i = 0; i < _rows; ++i) for (int j = 0; j < _cols; ++j) tabu_list[i][j] = 0; }
};

#endif // !_SRC_STRUCT_H_

