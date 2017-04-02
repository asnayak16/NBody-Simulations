// N-BODY SIMULATION Using GPUs
// Exact pair interactions with Leapfrog Integration
// Written by Ashwin Nayak, asnayak[at]ucsd.edu
// ----------------------------------------------------------
// Graphics Card : NVIDIA GTX 960M
// ----------------------------------------------------------
// References : 
// GPU GEMS 3 Documentation
// CUDA Sample Programs and Documentation
// ----------------------------------------------------------
#include <stdio.h>
#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
// CUDA Libraries
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
//#include "math_functions.h"
// Timer Libraries
#include <chrono>
//#include <ctime>

// GLOBAL PARAMETERS
#define PI 4*atan(1.)
#define G 4.498279e-3

// GLOBAL SIMULATION PARAMETERS
#define NPAR 10000

// GPU PARAMETERS
#define BLOCK_SIZE 512

struct Particle {
	float m, x, y, z, v_x, v_y, v_z;
};
// Input Data from File
void input_data(Particle *p) {
	std::ifstream input("plummer_init.dat");
	for (int i = 0; i < NPAR; i++) {
		input >> p[i].x
			>> p[i].y
			>> p[i].z
			>> p[i].v_x
			>> p[i].v_y
			>> p[i].v_z;
		p[i].m = 1.0e11/NPAR;
	}
	input.close();
}

// Leapfrog Integration in Device
__global__ void leapfrog(Particle *p, float dt) {
	// Each particle decided by thread index 
	int i = blockDim.x*blockIdx.x + threadIdx.x;
	
	if (i<NPAR) {
		//// ---------------- Position Verlet algorithm -------------------
		//// Update position 1/2
		p[i].x += 0.5 * dt * p[i].v_x;
		p[i].y += 0.5 * dt * p[i].v_y;
		p[i].z += 0.5 * dt * p[i].v_z;

		// Compute Acceleration
		float a_x = 0.; float a_y = 0.; float a_z = 0.;
			for (int j = 0; j < NPAR; j++) {
				if (i != j) {
					float r_x = p[i].x - p[j].x;
					float r_y = p[i].y - p[j].y;
					float r_z = p[i].z - p[j].z;
					float r_2 = r_x*r_x + r_y*r_y + r_z*r_z;
					float inv_r = rsqrtf(r_2); 
					a_x -=  p[j].m * r_x * inv_r * inv_r * inv_r;
					a_y -=  p[j].m * r_y * inv_r * inv_r * inv_r;
					a_z -=  p[j].m * r_z * inv_r * inv_r * inv_r;
				}
			}
		// Update Velocity
			p[i].v_x += dt * G * a_x;
			p[i].v_y += dt * G * a_y;
			p[i].v_z += dt * G * a_z;

		// Update Position 2/2
			p[i].x += 0.5 * dt * p[i].v_x;
			p[i].y += 0.5 * dt * p[i].v_y;
			p[i].z += 0.5 * dt * p[i].v_z;
		// --------------------------------------------------------------
	 } 
}

void output_data(Particle *p, int n) {
	std::string file_name = "pl_";

	file_name += std::to_string(n) + ".dat";
	std::ofstream output(file_name.c_str());
	//output  << p[i].m << "\t"
	for (int i = 0; i < NPAR; i++) {
		output << p[i].x << "\t"
			<< p[i].y << "\t"
			<< p[i].z << "\t"
			<< p[i].v_x << "\t"
			<< p[i].v_y << "\t"
			<< p[i].v_z << "\n";

	}
}

void PressEnterToContinue() {
	std::cout << "Press ENTER to continue... ";
	std::cin.ignore(std::numeric_limits <std::streamsize> ::max(), '\n');
}

// MAIN PROGRAM
int main() {	
	
	int nt, dsnap;
	float t, dt;
	cudaError_t cudaStatus;

	// Allocate Arrays
	int num_bytes = NPAR * sizeof(Particle);
	Particle *p = new Particle[NPAR];

	// Read Input Parameters from file 
	input_data(p);
	
	// Allocate Memory on Device
	cudaStatus = cudaSetDevice(0);
	Particle *dev_p = new Particle[NPAR];
	cudaStatus = cudaMalloc(&dev_p, num_bytes);
	if (cudaStatus != cudaSuccess) fprintf(stderr, "cudaMalloc failed!");

	int NBlocks = (NPAR + BLOCK_SIZE - 1) / BLOCK_SIZE;
	// Simulation Parameters 
	t = 0.;		dt = 0.005;
	nt = 3000;	dsnap = 12;
	
	// Start Timer
	auto start = std::chrono::system_clock::now();

	// Time Loop
	for (int it = 0; it < nt; it++) {
		t = t + dt;

		// Display Progress
		if (it%100==0)
		std::cout << "Iteration:" << it+1 << "\t Time: " << t << std::endl;
		
		// Copy memory to device
		cudaStatus = cudaMemcpy(dev_p, p, num_bytes, cudaMemcpyHostToDevice);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "cudaMemcpy failed!");
		
		// March all particles in time
		leapfrog<<<NBlocks,BLOCK_SIZE>>>(dev_p,dt);
		
		// Copy memory back to Host
		cudaStatus = cudaMemcpy(p, dev_p, num_bytes,cudaMemcpyDeviceToHost);
		if (cudaStatus != cudaSuccess) fprintf(stderr, "cudaMemcpy2 failed!");
		
		// Output every dsnap iterations
		if ((it+1)%dsnap == 0) output_data(p, it+1);
		
	}
	
	// Stop Timer
	auto finish = std::chrono::system_clock::now();
	
	// Display Timer info
	std::cout << "Time elapsed : " 
		<< std::chrono::duration_cast<std::chrono::seconds>(finish-start).count()/60.0
		<< " minutes\n";
	std::cout << "Avg time per iteration : "
		<< (std::chrono::duration_cast<std::chrono::seconds>(finish - start).count())/(float)nt
		<< " seconds\n";

	cudaFree(dev_p);
	delete[] p;

	PressEnterToContinue(); // Windows user's dilemma
	return 0;
}
