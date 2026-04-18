#include "hip/hip_runtime.h"
#include "hiprand/hiprand_kernel.h"


__global__ void wl_sampling(int* in, int* out, int M, int N,unsigned long long seed,int sampling)
{
    //if(blockIdx.x * blockDim.x + threadIdx.x>=M) return; //exit if exceeds the limits of matrix
    //if(blockIdx.y * blockDim.y + threadIdx.n>=N) return; 
    int idx = (blockIdx.x * blockDim.x + threadIdx.x) + (blockIdx.y * blockDim.y +
                threadIdx.y)*M;

    
    hiprandStateMRG32k3a_t state; //random number generator
    hiprand_init(seed,idx,0,&state); //initializes state and seed

    int holder[sampling]; //holds the number of energy changes
    for(int i =0;i<sampling;i++)
    {
        holder[i]=-1;
    }
    int elem=0; //number of elements in holder array
    int en = 0; //change in energy
    for(int i =0;i<sampling;i++)
    {
        bool binary = 0; //check for if value is already present in holder
        float r = hiprand_uniform(&state); //generates a randome number
        int x = int(r*M);
        r = hiprand_uniform(&state);
        int y = int(r*N);
        int check = x+y*M; //new index from random number generator
        for(int j =0;j<i;j++)
        {
            if(holder[j]==check)
            {
                binary=1; // check for presence of 
            }
        }
        if(binary==0)
        {
            holder[elem] = x + y*M; //adds element to list
            elem = elem + 1; //increments number of elements in list
            int c_idx = x + y*M; //current index 
        
            int n_u =  x + (y-1)*M; // defines the standard upper neighbor
            int n_d = (x) + ((y) +1)*M; // defines the standard lower neighbor
            int n_l = ((x)-1) + ((y))*M; // defines the standard left neighbor
            int n_r = ((x)+1) + ((y))*M; //defines the standard right neighbor 
            
            if(x == 0) // leftmost column pbc
            {
                n_l=(M-1) + ((y) -1)*M;
            }
            if(x == M-1) // rightmost column pbc
            {
                n_r=0 + ((y) -1)*M;
            }
            if(y == 0) // uppermost row pbc
            {
                n_u= (x) + (N-1)*M;
            }
            if(y == N-1) // lowermost row pbc
            {
                n_d= (x) + (0)*M;
            }

            // calculates if neighbor has been seen in holder array 
            int n_check = 1;
            for(int j =0;j<elem;j++) 
            {
                if(holder[j]==n_l)
                {
                    en+=-(in[c_idx])*-(in[n_l]);
                    n_check  =0;
                }
            }
            en+=(in[n_l]*-(in[c_idx]))*n_check;
            n_check  =1;
            for(int j =0;j<elem;j++)
            {
                if(holder[j]==n_r)
                {
                    en+=-(in[c_idx])*-(in[n_l]);
                    n_check  =0;
                }
            }
            en+=(in[n_r]*-(in[c_idx]))*n_check;
            n_check  =1;
            
            for(int j=0;j<elem;j++)
            {
                if(holder[j]==n_u)
                {
                    en+=-(in[c_idx])*-(in[n_u]);
                    n_check  =0;
                }
            }
            en+=(in[n_u]*-(in[c_idx]))*n_check;
            n_check  =1;
            for(int j=0;j<elem;j++)
            {
                if(holder[j]==n_d)
                {
                    en+=-(in[c_idx])*-(in[n_d]);
                    n_check  =0;
                }
            }
            en+=(in[n_d]*-(in[c_idx]))*n_check;
            // calculates if neighbor has been seen in holder array 
            
            
            
            printf("idx: %i current_en: %i neighbors: %i, %i, %i, %i  x: %i  y: %i cidx: %i",idx,en,n_u,n_d,n_l,n_r,x,y,c_idx);

        }
        /*
        int idx = (blockIdx.x * blockDim.x + threadIdx.x) + (blockIdx.y * blockDim.y +
                threadIdx.y)*M;
        
        int n_u =  (blockIdx.x * blockDim.x + threadIdx.x) + ( (blockIdx.y * blockDim.y +
                threadIdx.y) -1)*M;
        int n_d = (blockIdx.x * blockDim.x + threadIdx.x) + ( (blockIdx.y * blockDim.y +
                threadIdx.y) +1)*M;
        int n_l = ((blockIdx.x * blockDim.x + threadIdx.x)-1) + ( (blockIdx.y * blockDim.y +
                threadIdx.y))*M;
        int n_r = ((blockIdx.x * blockDim.x + threadIdx.x)+1) + ( (blockIdx.y * blockDim.y +
                threadIdx.y))*M;
        
        if(blockIdx.x * blockDim.x + threadIdx.x == 0)
        {
            n_l=(M-1) + ((blockIdx.y * blockDim.y +
                threadIdx.y) -1)*M;
        }
        if(blockIdx.x * blockDim.x + threadIdx.x == M-1)
        {
            n_r=0 + ((blockIdx.y * blockDim.y +
                threadIdx.y) -1)*M;
        }
        if(blockIdx.y * blockDim.y + threadIdx.y == 0)
        {
            n_u= (blockIdx.x * blockDim.x + threadIdx.x) + (N-1)*M
        }
        if(blockIdx.y * blockDim.y + threadIdx.y == N-1)
        {
            n_d= (blockIdx.x * blockDim.x + threadIdx.x) + (0)*M
        }
        */
        
    }

    
}

int main(int argc, char** argv)
{
    int nx = 4;//32;
    int ny = 4;//32;
    if(argc > 1) nx = atoi(argv[1]);
    if(argc > 2) ny = atoi(argv[2]);

    unsigned int M = 20;//4960;
    unsigned int N = 20;//4960;

    if(argc > 3) M = atoi(argv[3]);
    if(argc > 4) N = atoi(argv[4]);


    std::cout << "M: " << M << " N: " << N << std::endl;
    size_t size   = sizeof(int) * M * N;
    int*   matrix = (int*) malloc(size);

    int*   neighbor_change = (int*) malloc(size);
    int*   checked_location = (int*) malloc(size*2);

    for(int i = 0; i < M * N; i++)
    {
        //matrix[i] = rand() % 1002;
        matrix[i] = 1;
    }
    int *in, *out;

    hipMalloc(&in, size);
    hipMalloc(&out, size);
    hipMemset(in, 1, size);
    hipMemset(out, 1, size);
    hipMemcpy(in, matrix, size, hipMemcpyHostToDevice);
    hipDeviceSynchronize();
    hipDeviceProp_t props;
    hipGetDeviceProperties(&props, 0);


    dim3 grid(((M-1) / nx)+1, ((N-1) / ny)+1, 1);
    dim3 block(nx, ny, 1);  

    hipLaunchKernelGGL(wl_sampling, grid, block, 0, 0, in, out, M, N,1234ULL,1);
    check_hip_error();
    return 0;


    hipFree(in);
    hipFree(out);

    free(matrix);

    return 0;
}
