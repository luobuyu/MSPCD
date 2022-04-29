#ifndef _SRC_UTILS_H_
#define _SRC_UTILS_H_

#include <iostream>
#include <cstdio>
#include <string>
#include <sstream>
#include <fstream>
#include <ctime>
#include <iomanip>
#include <random>
#include <chrono>
using namespace std;
namespace utils {

	using namespace std;

	class Date {
	public:
		// 返回表示日期格式的字符串，年月日
		static string to_short_str() {
			time_t now = time(nullptr);
			char mbstr[100];
			strftime(mbstr, sizeof(mbstr), "%y%m%d", localtime(&now));
			return mbstr;
		}
		// 返回表示日期格式的字符串，年月日时分秒
		static string to_long_str() {
			time_t now = time(nullptr);
			char mbstr[100];
			strftime(mbstr, sizeof(mbstr), "%y%m%d%H%M%S", localtime(&now));
			return mbstr;
		}
		// 返回日期的格式化字符串，年/月/日 时:分:秒
		static string to_format_str()
		{
			time_t now = time(nullptr);
			char mbstr[100];
			strftime(mbstr, sizeof(mbstr), "%y/%m/%d %H:%M:%S", localtime(&now));
			return mbstr;
		}
	};

	class Timer {
	private:
		using Clock = chrono::steady_clock;
		using TimePoint = Clock::time_point;
		TimePoint startTime, endTime;
	public:
		Timer(long long secTimeout)
		{
			setTimer(secTimeout);
		}
		void setTimer(long long secTimeout)
		{
			startTime = Clock::now();
			endTime = startTime + chrono::seconds(secTimeout);
		}
		bool isTimeout()
		{
			return endTime < Clock::now();
		}
		double getDuration()
		{
			return std::chrono::duration<double>(Clock::now() - startTime).count();
		}
	};
	
	static void split_filename(const string& str, string& dir, string& file, string& id) {
		size_t found1 = str.find_last_of("/\\");
		size_t found2 = str.find_last_of("_");
		size_t found3 = str.find_last_of(".");
		// substr，起点，个数
		dir = str.substr(0, found1 + 1);
		file = str.substr(found1 + 1, found3 - found1 - 1);
		id = found2 == string::npos ? string() : str.substr(found2, found3 - found2);
	}
}


namespace fast_io
{
    inline char nextChar()
    {
        static char buf[1000000], * p1 = buf, * p2 = buf;
        return p1 == p2 && (p2 = (p1 = buf) + fread(buf, 1, 1000000, stdin), p1 == p2) ? EOF : *p1++;
    }
#define getch getchar
    template <class T>
    inline void read(T& x)
    {
        T flag = 1;
        x = 0;
        char ch = getch();
        while (ch < '0' || ch > '9')
        {
            if (ch == '-')
                flag = -1;
            ch = getch();
        }
        while (ch >= '0' && ch <= '9')
        {
            x = (x << 3) + (x << 1) + (ch ^ 48), ch = getch();
        }
        x *= flag;
    }

    template <class T, class... _T>
    inline void read(T& x, _T &...y)
    {
        return read(x), read(y...);
    }

    template <class T>
    inline void print128(T x)
    {
        if (x < 0)
            putchar('-'), x = -x;
        if (x > 9)
            print128(x / 10);
        putchar(x % 10 + '0');
    }

	enum class LineType
	{
		HEAD, LINE
	};

	template <class T>
	inline void writeStream(T x)
	{
		cerr << x << endl;
	}

	template <class T, class... _T>
	inline void writeStream(T x, _T ...y)
	{
		cerr << x << ",";
		writeStream(y...);
	}

	// 写入日志文件，参数 路径，表头还是行数据，具体内容
	// 文件为空才会创建表头，否则的话
	template <class T, class... _T>
	inline void writeLog(string logPath, T x, _T ...y)
	{
		freopen(logPath.c_str(), "a", stderr);
		ofstream log_file(logPath, ios::app);
		log_file.seekp(0, ios::end);
		if (is_same<T, LineType>::value && x == LineType::HEAD)
		{
			if (log_file.tellp() <= 0)
			{
				writeStream(y...);
			}
		}
		else if (is_same<T, LineType>::value && x == LineType::LINE)
		{
			writeStream(y...);
		}
	}
	template <class T, class... _T>
	inline void writeLog(T x, _T ...y)
	{
		if (is_same<T, LineType>::value)
		{
			writeStream(y...);
		}
	}

} // namespace FAST_IO

namespace Debug {
#define DEBUG
#ifdef DEBUG
#define eprintf(x) fprintf(stderr, x)
#define eputs(str) fputs(str, stderr), fputc('\n', stderr)
#define var(x) ""#x" = " << x
#define watch(...) trace(#__VA_ARGS__, __VA_ARGS__)
#else
#define eprintf(...)
#define eputs(...)
#define watch(...)
#endif

	template <class printable>
	void trace(const char* name, printable&& value)
	{
		cerr << name << " = " << value << endl;
	}
	template <class printable, class ...args>
	void trace(const char* names, printable&& value, args &&...list)
	{
		const char* separate = strchr(names + 1, ',');
		cerr.write(names, separate - names) << " = " << value << ',';
		trace(separate + 1, list...);
	}
}

namespace utils_visualize_drawer {

	using namespace std;

	class Random {
	public:
		using Generator = std::mt19937;


		Random(int seed) : rgen(seed) {}
		Random() : rgen(generateSeed()) {}


		static int generateSeed() {
			return static_cast<int>(std::time(nullptr) + std::clock());
		}

		Generator::result_type operator()() { return rgen(); }

		// pick with probability of (numerator / denominator).
		bool isPicked(unsigned numerator, unsigned denominator) {
			return ((rgen() % denominator) < numerator);
		}

		// pick from [min, max).
		int pick(int min, int max) {
			return ((rgen() % (max - min)) + min);
		}
		// pick from [0, max).
		int pick(int max) {
			return (rgen() % max);
		}

	protected:
		Generator rgen;
	};

	struct RandColor {
		static constexpr auto ColorCodeChar = "0123456789ABCDEF";
		static constexpr int ColorCodeBase = 16;
		static constexpr int ColorCodeLen = 6;

		void next() {
			for (int i = 0; i < ColorCodeLen; ++i) {
				int c = r.pick(ColorCodeBase);
				bcolor[i] = ColorCodeChar[c];
				fcolor[i] = ColorCodeChar[(c > (ColorCodeBase / 2)) ? 0 : (ColorCodeBase - 1)]; // (c + ColorCodeBase / 2) % ColorCodeBase
			}
		}

		char fcolor[ColorCodeLen + 1] = { 0 }; // front color.
		char bcolor[ColorCodeLen + 1] = { 0 }; // background color.
		Random r;
	};

	struct Drawer {
		static constexpr double W = 400;
		static constexpr double H = 300;

		Drawer(string path, double width, double height) : ofs(path), wx(W / width), hx(H / height) { begin(); }
		~Drawer() { end(); }

		void begin() {
			ofs << "<!DOCTYPE html>" << endl
				<< "<html>" << endl
				<< "  <head>" << endl
				<< "    <meta charset='utf-8'>" << endl
				<< "    <title>SmartMPW Visualization</title>" << endl
				<< "  </head>" << endl
				<< "  <body>" << endl // style='text-align:center;'
				<< "    <svg width='" << W << "' height='" << H << "' viewBox='-50 -50 " << W + 100 << " " << H + 100 << "'>" << endl;
		}
		void end() {
			ofs << "    </svg>" << endl
				<< "  </body>" << endl
				<< "</html>" << endl;
		}

		void rect(double x, double y, double w, double h, bool d, const string& label, const string& fcolor, const string& bcolor) {
			if (d) { swap(w, h); }
			x *= wx; y *= hx; w *= wx; h *= hx;
			ofs << "      <rect x='" << x << "' y='" << y << "' width='" << w << "' height='" << h << "' style='fill:#" << bcolor << "; stroke:black; stroke-width:2'/>" << endl
				<< "      <text x='" << x + w / 2 << "' y='" << y + h / 2 << "' text-anchor='middle' alignment-baseline='middle' style='fill:#" << fcolor << "'>" << label << "</text>" << endl << endl;
		}
		void rect(double x, double y, double w, double h, bool d = false, const string& label = "") {
			rc.next();
			rect(x, y, w, h, d, label, rc.fcolor, rc.bcolor);
		}

		void wire(double x, double y, double w, double h, const string& label = "") {
			rc.next();
			x *= wx; y *= hx; w *= wx; h *= hx;
			ofs << "      <rect x='" << x << "' y='" << y << "' width='" << w << "' height='" << h << "' style='fill:none" << "; stroke:#" << rc.bcolor << "; stroke-width:2; stroke-dasharray:12, 4'/>" << endl
				<< "      <text x='" << x + w / 2 << "' y='" << y + h / 2 << "' text-anchor='middle' alignment-baseline='middle' style='fill:#" << rc.fcolor << "'>" << label << "</text>" << endl << endl;
		}

		void line(double x1, double y1, double x2, double y2, int layer) {
			static vector<int> cutWidth({ 8, 8, 6 });
			static vector<string> cutColor({ "red", "blue", "orange" });
			x1 *= wx; y1 *= hx; x2 *= wx; y2 *= hx;
			ofs << "      <line x1='" << x1 << "' y1='" << y1 << "' x2='" << x2 << "' y2='" << y2 << "' stroke-dasharray='12, 4' stroke='" << cutColor[layer] << "' stroke-width='" << cutWidth[layer] << "'/>" << endl << endl;
		}

		void circle(double x, double y, double r = 2) {
			x *= wx; y *= hx;
			ofs << "      <circle cx='" << x << "' cy='" << y << "' r='" << r << "' style='fill-opacity:0; stroke:#000000; stroke-width:2'/>" << endl << endl;
		}

		void polygon(const string& polygon_str, const string& label, const string& fcolor, const string& bcolor) {
			ofs << "      <polygon points='" << polygon_str << "' style='fill:#" << bcolor << "; stroke:black; stroke-width:1'/>" << endl << endl;
		}
		void polygon(const string& polygon_str, const string& label = "") {
			rc.next();
			polygon(polygon_str, label, rc.fcolor, rc.bcolor);
		}


		double wx;
		double hx;
		ofstream ofs;
		RandColor rc;
	};
}


#endif // !_SRC_UTILS_H_

