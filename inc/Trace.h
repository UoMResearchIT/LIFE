
#pragma once

#include <string>
#include <iostream>
#include <fstream>
//#include <ifstream>
#include <functional>
#include <sstream>

template <typename TT>
std::string number_to_string(const TT &value)
{
	std::stringstream ss;
	ss.precision(17);
	ss << value;
	return ss.str();
}

class Trace
{
public:

	Trace(bool enabled=false) : enabled(enabled) {
		std::cout << "enabled = " << enabled << std::endl;
		if (enabled && fileExists(out_filename)) {
			// Truncate any existing file
			std::ofstream os;
			os.open(out_filename, std::ofstream::out);
		}
	}

	void addTracer(const std::string &name, std::function<std::string(void)> f) {
		if (!enabled) {return;};
		for (auto &tracer : tracers) {
			if (tracer.first == name) {
				tracer.second = f;
				return;
			}
		}
		tracers.push_back(std::make_pair(name,f));
		header();
	};

	void removeTracer(const std::string &name) {
		if (!enabled) {return;};
		for (auto &tracer : tracers) {
			if (tracer.first == name) {
				tracer.second = [](){return "-";};
				return;
			}
		}
		assert(0);
	};

	
	void header() {
		if (!enabled) {return;};
		std::ofstream os;
		os.open(out_filename, std::ofstream::out | std::ofstream::app);

		os << "Tracepoint\t";

		for (auto &tracer : tracers) {
			os << tracer.first << "\t";
		}
		os << std::endl;
	}
	void trace(const std::string &traceName) {
		if (!enabled) {return;};
		std::ofstream os;
		os.open(out_filename, std::ofstream::out | std::ofstream::app);

		os << traceName << "\t";
		for (auto &tracer : tracers) {
			os << tracer.second() << "\t";
		}
		os << std::endl;
	}

	private:

	bool fileExists(const std::string &filename) {
		std::ifstream f(filename.c_str());
		std::cout << "Trace file already exists? " << f.good() << std::endl;
		return f.good();		
	}

	std::vector<std::pair<std::string,std::function<std::string(void)>>> tracers;
	std::string out_filename = "trace.txt";
	bool enabled;
};


extern Trace trace;

class Tracer
{
	public:

	Tracer(const std::string &name, std::function<std::string(void)> f) {
		trace.addTracer(name,f);
		this->name = name;
	}

	~Tracer() {
		trace.removeTracer(name);
	}

	private:

	std::string name;
};
