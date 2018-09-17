#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include<stdbool.h>
#include <math.h>
#include <omp.h>
#include<time.h>
float distance2(const float *v1, const float *v2, const int d) {
  float dist = 0.0;
  for (int i=0; i<d; i++) {
    float diff = v1[i] - v2[i];
    dist += diff * diff;
  }
  return dist;
}
double cpu_time(void)
{
  double value;
  value=(double) clock()/ (double) CLOCKS_PER_SEC;
  return value;
}
float* search(float *site2,float *centroids,int k,int d)
{
	//int *closest;
	printf("NEarest\n");
	double min,dist;
	int cent,m=0;
	int n=sizeof(site2);
	for(int i=0;i<k;i++)
	{
		min=99999.00;
		for(int j=0;j<n;j++)
		{
			dist=0.0;
			for(int l=0;l<d;l++)
			{
				dist+=pow(site2[j*d+l]-centroids[i*d+l],2);
			}
			dist=sqrt(dist);
			if(dist<min)
			{
				min=dist;
				cent=j;
			}
		}
		for(int p=0;p<d;p++)
		{
			centroids[m]=site2[cent*d+p];
			m++;
		}
	}
	return centroids;
}
// Creates an array of random floats. Each number has a value from 0 - 1
float* create_rand_nums(const int num_elements) {
  float *rand_nums = (float *)malloc(sizeof(float) * num_elements);
  assert(rand_nums != NULL);
  for (int i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() / (float)(RAND_MAX) );
  }
  return rand_nums;
}
float distance3(const float *v1, const float *v2,const int d)
{
  float dist=0.0;
  for(int i=0;i<d;i++)
  {
    dist+=abs(v1[i]-v2[i]);
  }
  return dist;
}
// Distance**2 between d-vectors pointed to by v1, v2.


// Assign a site to the correct cluster by computing its distances to
// each cluster centroid.
int assign_site(const float* site, float* centroids,
    const int k, const int d) {
  int best_cluster = 0;
  float best_dist = distance2(site, centroids, d);
  float* centroid = centroids + d;
  for (int c = 1; c < k; c++, centroid += d) {
    float dist = distance2(site, centroid, d);
    if (dist < best_dist) {
      best_cluster = c;
      best_dist = dist;
    }
  }
  return best_cluster;
}


// Add a site (vector) into a sum of sites (vector).
void add_site(const float * site, float * sum, const int d) {
  for (int i=0; i<d; i++) {
    sum[i] += site[i];
  }
}

// Print the centroids one per line.
void print_centroids(float * centroids, const int k, const int d) {
  float *p = centroids;
  printf("Centroids:\n");
  for (int i = 0; i<k; i++) {
    for (int j = 0; j<d; j++, p++) {
      printf("%f ", *p);
    }
    printf("\n");
  }
}

int main(int argc, char** argv) {
  if (argc != 4) {
    fprintf(stderr,
      "Usage: kmeans num_sites_per_proc num_means num_dimensions\n");
    exit(1);
  }
  double start=cpu_time();
  int m=0;
  // Get stuff from command line:
    // number of sites per processor.
    // number of processors comes from mpirun command line.  -n
  int sites_per_proc = atoi(argv[1]);
  int k = atoi(argv[2]);  // number of clusters.
  int d = atoi(argv[3]);  // dimension of data.
  // Seed the random number generator to get different results each time
  //  srand(time(NULL));
  // No, we'd like the same results.
  srand(31359);
  int mat[sites_per_proc][k];
  /*for(int i=0;i<sites_per_proc;i++)
    for(int j=0;j<k;j++)
      mat[i][j]=0;*/
  // Initial MPI and find process rank and number of processes.
  MPI_Init(NULL, NULL);
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  //
  // Data structures in all processes.
  //
  // The sites assigned to this process.
  float* sites;  
  assert(sites = malloc(sites_per_proc * d * sizeof(float)));
  // The sum of sites assigned to each cluster by this process.
  // k vectors of d elements.
  float* sums;
  assert(sums = malloc(k * d * sizeof(float)));
  // The number of sites assigned to each cluster by this process. k integers.
  int* counts;
  assert(counts = malloc(k * sizeof(int)));
  // The current centroids against which sites are being compared.
  // These are shipped to the process by the root process.
  float* centroids;
  assert(centroids = malloc(k * d * sizeof(float)));
    float* centroids2;
  assert(centroids2 = malloc(k * d * sizeof(float)));
  // The cluster assignments for each site.
  int* labels;
  assert(labels = malloc(sites_per_proc * sizeof(int)));
  
  //
  // Data structures maintained only in root process.
  //
  // All the sites for all the processes.
  // site_per_proc * nprocs vectors of d floats.
  float* all_sites = NULL;
  // Sum of sites assigned to each cluster by all processes.
  float* grand_sums = NULL;
  // Number of sites assigned to each cluster by all processes.
  int* grand_counts = NULL;
  // Result of program: a cluster label for each site.
  int* all_labels;
  if (rank == 0) {
    all_sites = create_rand_nums(d * sites_per_proc * nprocs);
     if ((rank == 0) && 1) {
    float* site = all_sites; 
    for (int i = 0;
   i < nprocs * sites_per_proc;
   i++, site += d) {
      for (int j = 0; j < d; j++) printf("%f ", site[j]);
      printf("\n");
    }
  }
    // Take the first k sites as the initial cluster centroids.
    for (int i = 0; i < k * d; i++) {
      centroids[i] = all_sites[i]; 
    }
    print_centroids(centroids, k, d);
    assert(grand_sums = malloc(k * d * sizeof(float)));
    assert(grand_counts = malloc(k * sizeof(int)));
    assert(all_labels = malloc(nprocs * sites_per_proc * sizeof(int)));
  }

  // Root sends each process its share of sites.
  MPI_Scatter(all_sites,d*sites_per_proc, MPI_FLOAT, sites,
              d*sites_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

  
  float norm = 1.0;  // Will tell us if centroids have moved.
  int itr=0;
  while (norm>0.01){ // While tthey've moved...

    // Broadcast the current cluster centroids to all processes.
    MPI_Bcast(centroids, k*d, MPI_FLOAT,0, MPI_COMM_WORLD);

    // Each process reinitializes its cluster accumulators.
    for (int i = 0; i < k*d; i++) sums[i] = 0.0;
    for (int i = 0; i < k; i++) counts[i] = 0;

    // Find the closest centroid to each site and assign to cluster.
    float* site = sites;
    for (int i = 0; i < sites_per_proc; i++, site += d) {
      int cluster = assign_site(site, centroids, k, d);
      // Record the assignment of the site to the cluster.
      counts[cluster]++;
      mat[i][cluster]=1;
      for(int p=0;p<k;p++)
      {
        if(p!=cluster)
          mat[i][p]=0;
      }
      add_site(site, &sums[cluster*d], d);
    }

    // Gather and sum at root all cluster sums for individual processes.
    MPI_Reduce(sums, grand_sums, k * d, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(counts, grand_counts, k, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    float *site2=sites;
    if (rank == 0) {

      // Root process computes new centroids by dividing sums per cluster
      // by count per cluster.
           for (int i = 0; i<k; i++) {
  for (int j = 0; j<d; j++) {
    int dij = d*i + j;
    grand_sums[dij] /= grand_counts[i];
  }
      }
      for(int i=0;i<k;i++)
      {int ctr=0;

        double dist,min2=99999.0;
        int cent;
        for(int j=0;j<sites_per_proc;j++)
        {
          dist=0.0;
          for(int t=0;t<sites_per_proc;t++)
          {
            if(mat[j][i]==1 && mat[t][i]==1)
            {
                dist+=distance2(site2+j*d,site2+t*d,d);
              
            }
          }
          if(dist<min2)
          {
          	if(cent!=j)
          	{
            	min2=dist;
            	cent=j;
            }
          }
        }
        	for(int x=0;x<d;x++)
        	{
        		centroids2[i*d+x]=sites[cent*d+x];
      
        	}
        
       
      // Have the centroids changed much?
      // Copy new centroids from grand_sums into centroids.
      
    }
    //print_centroids(centroids,k,d);
    norm = distance2(grand_sums, centroids, d*k);
      double change=0.0;
      for(int y=0;y<k*d;y++)
      {
          change+=(pow(centroids2[y]-centroids[y],2));
      }
      change=sqrt(change);
      printf("%f",change);
      if(norm>0.01)
        for(int u=0;u<k*d;u++)
        centroids[u]=grand_sums[u];
      printf("norm: %f\n",norm);
    
  		print_centroids(centroids,k,d);
    //  print_centroids(centroids,k,d);
    }
    // Broadcast the norm.  All processes will use this in the loop test.
    MPI_Bcast(&norm, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
  }
  if(rank==0)
   {centroids=search(sites,centroids,k,d);
	print_centroids(centroids,k,d);}
  // Now centroids are fixed, so compute a final label for each site.
  float* site = sites;
  for (int i = 0; i < sites_per_proc; i++, site += d) {
    labels[i] = assign_site(site, centroids, k, d);
  }
 
  // Gather all labels into root process.
  MPI_Gather(labels, sites_per_proc, MPI_INT,
       all_labels, sites_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

  // Root can print out all sites and labels.
  if ((rank == 0) && 1) {
    float* site = all_sites; 
    for (int i = 0;
   i < nprocs * sites_per_proc;
   i++, site += d) {
      for (int j = 0; j < d; j++) printf("%f ", site[j]);
      printf("%4d\n", all_labels[i]);
    }
  }
  if(rank==0)
  {
  double stop=cpu_time();  
  printf("Time elapsed: %f", stop-start);  
  }
  MPI_Finalize();

}