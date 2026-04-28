#include "hip/hip_runtime.h"
#include "hiprand/hiprand_kernel.h"


__global__ void
wl_sampling(int* in, int* out, int M, int N,unsigned long long seed,int sampling,
      int current_en ,int* hist;double* entropy_est)
{
    //if(blockIdx.x * blockDim.x + threadIdx.x>=M) return; //exit if exceeds the limits of matrix
    //if(blockIdx.y * blockDim.y + threadIdx.n>=N) return; //though this is deprecated now from algorithim changes

    int idx = (blockIdx.x * blockDim.x + threadIdx.x) + (blockIdx.y * blockDim.y + threadIdx.y)*M;

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
            
            
            if random() < exp(entropy[current_energy] - entropy[proposed_energy]):
                # If accepted, update the energy and the system:
                current_energy = proposed_energy
                system.accept_proposed_configuration()
            else:
                # If rejected
                system.reject_proposed_configuration()

            H[current_energy] += 1
            entropy[current_energy] += f

            printf("idx: %i current_en: %i neighbors: %i, %i, %i, %i  x: %i  y: %i cidx: %i",idx,en,n_u,n_d,n_l,n_r,x,y,c_idx);
            holder[elem] = x + y*M; //adds element to list
            elem = elem + 1; //increments number of elements in list
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
    
    if(idx==0)
    {
          
    }

    //int iidx = (blockIdx.x * blockDim.x + threadIdx.x) * N + (blockIdx.y * blockDim.y +
    //           threadIdx.y);
    //int oidx =
    //    (blockIdx.y * blockDim.y + threadIdx.y) * M + (blockIdx.x * blockDim.x + threadIdx.x);

    //out[oidx] = in[iidx];
}

__global__ void
energy_sum(int* in_matrix, int* en_out, int M, int N)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    int global_size = gridDim.x * blockDim.x;
    double local_sum = 0;

    if(idx>M*N-1) return;

    en_out[idx]=0;

    for (int i = idx; i < size; i += global_size) {

        if(i>M*N-1) continue;

        int n_u= idx - M; // defines the standard upper neighbor
        int n_d= idx + M; // defines the standard lower neighbor
        int n_l= idx - 1; // defines the standard left neighbor
        int n_r= idx + 1; //defines the standard right neighbor     

        if(n_u < 0) // leftmost column pbc
        {
            n_u = (M*(N-1))+idx;
        }
        if(n_d >= M*N) // leftmost column pbc
        {
            n_d = idx-(M*(N-1));
        }
        if(n_l%M==(M-1)) // rightmost column pbc
        {
            n_l=n_l+M;
        }
        if(n_r%M == 0) // uppermost row pbc
        {
            n_r=n_r-M;
        }
        
        en+=(in[n_u]*(in[c_idx]));
        en+=(in[n_d]*(in[c_idx]));
        en+=(in[n_l]*(in[c_idx]));
        en+=(in[n_r]*(in[c_idx]));
        en_out[idx] += input[i];
    }

    //atomicAdd(output,local_sum);        
    
}

int energy_rocprim(int &en_matrix, int M, int N)
{
    size_t input_size = N*M*sizeof(int);    // e.g., 8
    //short * input;        // e.g., [1, 2, 3, 4, 5, 6, 7, 8]
    int * output;         // empty array of 1 element

    size_t temporary_storage_size_bytes;
    void * temporary_storage_ptr = nullptr;
    // Get required size of the temporary storage
    rocprim::reduce(
        temporary_storage_ptr, temporary_storage_size_bytes,
        en_matrix, output, input_size, rocprim::plus<int>()
    );

    // allocate temporary storage
    hipMalloc(&temporary_storage_ptr, temporary_storage_size_bytes);

    // perform reduce
    rocprim::reduce(
        temporary_storage_ptr, temporary_storage_size_bytes,
        en_matrix, output, input_size, rocprim::plus<int>()
    );

    hipDeviceSynchronize(temporary_storage_ptr);

    return output[0];
}

int
main(int argc, char** argv)
{
    int samples = 10;
    int en = 0;
    int d_en = 0;

    int nx = 4;//32;
    int ny = 4;//32;
    if(argc > 1) nx = atoi(argv[1]);
    if(argc > 2) ny = atoi(argv[2]);

    unsigned int M = 20;//4960;
    unsigned int N = 20;//4960;

    if(argc > 3) M = atoi(argv[3]);
    if(argc > 4) N = atoi(argv[4]);


    std::cout << "M: " << M << " N: " << N << std::endl;
    size_t size = M * N;

    dim3 grid(((M-1) / nx)+1, ((N-1) / ny)+1, 1); //default (128,128,1)
    dim3 block(nx, ny, 1);  //default (32,32,1) // transpose_a
    // dim3 grid(M/64, N/64, 1); dim3 block(64, 8, 1); // transpose_e


    int*   matrix = (int*) malloc(size*sizeof(int));
    int*   en_matrix = (int*) malloc(size*sizeof(int));

    double f_factor = 1;
    int* histogram = (int*) malloc(10*size*sizeof(int));
    double* estimation = (double*) malloc(10*size*sizeof(double));

    int*   neighbors_to_change = (int*) malloc(samples*sizeof(int));
    for(int i = 0; i < samples; i++)
    {
         neighbors_to_change[i] = -1;
    }

    for(int i = 0; i < size*10; i++)
    {
        histogram[i] = 0;
        estimation[i] = 0;
    }
    for(int i = 0; i < M * N; i++)
    {
        matrix[i] = 1;
    }
    int *d_in_matrix, *d_out_matrix, *d_en_matrix;
    

    std::chrono::high_resolution_clock::time_point t1, t2;

    hipMalloc(&d_in_matrix, size*sizeof(int));
    //hipMalloc(&d_out_matrix, size);

    hipMalloc(&d_en_matrix, size*sizeof(int)); //energy matrix


    double* d_estimation; //device wl estimations
    int* d_histogram;
    double* d_f_factor;
    hipMalloc(&d_estimation, size*10*sizeof(double));
    hipMalloc(&d_histograms, size*10*sizeof(int));
    //hipMalloc(&d_f_factor, sizeof(double));
    
    int* d_neighbors_to_change;
    hipMalloc(&d_neighbors_to_change, 10*sizeof(int));
    hipMemset(d_neighbors_to_change,-1, 10*sizeof(int));

    hipMemset(d_in_matrix, 0, size*sizeof(int));
    //hipMemset(d_out_matrix, 0, size);
    hipMemset(d_en_matrix, 0, size*sizeof(int));

    hipMemset(d_estimation, 0,size*10*sizeof(double));
    hipMemset(d_histograms, 0,size*10*sizeof(int));

    check_hip_error();
    hipMemcpy(d_in_matrix, matrix, size, hipMemcpyHostToDevice);
    hipDeviceSynchronize();
    check_hip_error();
    hipDeviceProp_t props;
    hipGetDeviceProperties(&props, 0);

    dim3 grid(((M-1) / nx)+1, ((N-1) / ny)+1, 1); //default (128,128,1)
    dim3 block(nx, ny, 1);  //default (32,32,1) // transpose_a
    // dim3 grid(M/64, N/64, 1); dim3 block(64, 8, 1); // transpose_e

    hipLaunchKernelGGL(energy_sum, grid, block, 0, 0, d_in_matrix, d_en_matrix, M, N);

    hipMemcpy(en_matrix, d_en_matrix, size, hipMemcpyDeviceToHost);
    
    energy_rocprim(en_matrix,M,N);
    // warmup
    hipLaunchKernelGGL(wl_sampling, grid, block, 0, 0, in, out, M, N,1234ULL,samples,en,d_histograms,d_estimation,d_en);

    


    check_hip_error();
    return 0;
}
