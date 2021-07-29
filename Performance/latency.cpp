// ConsoleApplication2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <array>
#include<chrono>
//#include "ConsoleApplication2.h"
#include <vector>
#include<numeric>

long long ReadLatency()
{
	auto t1 = std::chrono::high_resolution_clock::now();

	//int size = 1024 * 1024 * 1024;//1 GB
	int size = 1024 * 1024;//1 MB
	int sum = 0;
	char* arr = (char*)malloc(size);
	for (int i = 0; i < size; i++)
		sum += arr[i];//access the memory and do some operation
	auto t2 = std::chrono::high_resolution_clock::now();
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();

	std::cout << "Duration in microseconds: " << duration << '\n';
	//std::cout << "Print the sum: " << sum << '\n';
	return duration;
}

//run the code with optimizations turned off
int main()
{
	int runs = 10;
	std::vector<long long> timings(runs);
	for (size_t i = 0; i < runs; i++)
	{
		timings.push_back(ReadLatency());
	}
	int average = 0;
	average = std::accumulate(timings.begin(), timings.end(), 0) / runs;
	std::cout << "Average duration in microseconds: " << average << '\n';
}
