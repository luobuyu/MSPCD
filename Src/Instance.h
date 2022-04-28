#ifndef _SRC_INSTANCE_H_
#define _SRC_INSTANCE_H_

#include <iostream>
#include <list>
#include <fstream>
#include <string>
#include "Utils.h"
#include "Config.h"
#include "Struct.h"
using namespace std;

class Environment {

private:
	string _ins_path;
	string _ins_dir;
	string _ins_name;
	string _ins_id;

public:
	//Environment(const string& ins_str) : _ins_path(ins_str) { utils::split_filename(_ins_path, _ins_dir, _ins_name, _ins_id); }
	void set_instance_name(string name) { _ins_name = name; }
public:
	const string& instance_name() const { return _ins_name; }
	string instance_path() const { return instance_dir() + _ins_name + file_suffix(); }
	string solution_path() const { return solution_dir() + _ins_name + file_suffix(); }
	string solution_path_with_time() const { return solution_dir() + _ins_name + "." + utils::Date::to_long_str() + file_suffix(); }
	string ins_html_path() const { return instance_dir() + _ins_name + ".html"; }
	string sol_html_path() const { return solution_dir() + _ins_name + ".html"; }
	string sol_html_path_with_time() const { return solution_dir() + _ins_name + "." + utils::Date::to_long_str() + ".html"; }
	string log_path() const { return solution_dir() + "log.csv"; }
private:
	static string file_suffix() { return ".td"; }
	static string instance_dir() { return "../Deploy/Instance/Bench/"; }
	static string solution_dir() { return "../Deploy/Solution/Bench/"; }
};

// 读取样例，并将样例处理后存起来
class Instance {

private:
	const Environment& _env; // 样例的配置路径
	int _task_num;
	int _vir_task_num;		// 等于 _task_num + 2
	vector<Task> _tasks;
	vector<vector<Edge>> _js;
	vector<vector<Edge>> _jp;
public:
	Instance(const Environment& env) : _env(env) { read_instance(); }

	/*********getter and setter************/
	vector<Task>& get_tasks() { return _tasks; }
	vector<vector<Edge>>& get_js() { return _js; }
	vector<vector<Edge>>& get_jp() { return _jp; }
	int get_task_num() { return _task_num; }
	int get_vir_task_num() { return _vir_task_num; }

private:
	void read_instance() {
		cin >> _task_num;
		_vir_task_num = _task_num + 2;
		int start = 0, end = _task_num + 1;
		// 0, 1 - _task_num, _task_num + 1
		_tasks.resize(_vir_task_num);
		_js.resize(_vir_task_num);
		_jp.resize(_vir_task_num);
		for (int i = 1, u, length, s_num; i <= _task_num; ++i)
		{
			cin >> u >> length >> s_num;
			_tasks[u] = Task(length);
			//_tasks[u].length = length;
			for (int j = 1, v, w; j <= s_num; ++j)
			{
				cin >> v >> w;
				_js[u].emplace_back(u, v, w);
				_jp[v].emplace_back(v, u, w);
			}
		}
		_tasks[start] = _tasks[end] = Task(0);
		//_tasks[start].length = _tasks[end].length = 0;
		for (int i = 1; i <= _task_num; ++i)
		{
			if (_js[i].size() == 0)
			{
				// 没有后缀，说明需要连接 end
				_js[i].emplace_back(i, end, 0);
				_jp[end].emplace_back(end, i, 0);
			}
			if (_jp[i].size() == 0)
			{
				// 没有前缀，说明需要连接 start
				_js[start].emplace_back(start, i, 0);
				_jp[i].emplace_back(i, start, 0);
			}
		}
		//visual();
	}

	void visual()
	{
		cerr << _vir_task_num << endl;
		for (int i = 0; i < _vir_task_num; ++i)
		{
			for (auto& item : _js[i])
			{
				cerr << item.u << " " << item.v << " " << item.w << endl;
			}
		}
	}
};

#endif // !_SRC_INSTANCE_H_

