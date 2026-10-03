// Exact float64 k-nearest neighbors, with certified nonnegative-distance pruning.
// Build: c++ -O3 -std=c++17 -shared -fPIC nearest.cpp -o nearest.so
// No fast-math: stable scientific floating-point/tie semantics are retained.
#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

extern "C" void metric_scores(int n, int nt, int q, int d, const int* order,
    const double* mol, const double* hd, const int* permtest, const int* permtrain,
    const double* yt, const double* actual, double lambda, double* result) {
    const int ks[3]={10,25,50};
    std::fill(result,result+q*3,0.0);
    for(int i=0;i<n;++i) {
        std::vector<std::pair<double,int>> heap; heap.reserve(51);
        for(int r=0;r<nt;++r) {
            int j=order[i*nt+r]; double m=mol[i*nt+j];
            if(heap.size()==50 && m>heap.front().first) break;
            double full=m+lambda*hd[permtest[i]*nt+permtrain[j]];
            std::pair<double,int> v(full,j);
            if(heap.size()<50) { heap.push_back(v); std::push_heap(heap.begin(),heap.end()); }
            else if(v<heap.front()) {
                std::pop_heap(heap.begin(),heap.end()); heap.back()=v;
                std::push_heap(heap.begin(),heap.end());
            }
        }
        std::sort(heap.begin(),heap.end());
        for(int t=0;t<q;++t) {
            std::vector<double> sum(d,0.0); int z=0;
            double anorm=0;
            for(int a=0;a<d;++a) anorm+=actual[(i*q+t)*d+a]*actual[(i*q+t)*d+a];
            for(int r=0;r<50;++r) {
                int j=heap[r].second;
                for(int a=0;a<d;++a) sum[a]+=yt[(j*q+t)*d+a];
                if(r+1==ks[z]) {
                    double norm=0; for(double v:sum) norm+=v*v;
                    norm=std::sqrt(norm); if(norm==0) norm=1;
                    double error=0;
                    for(int a=0;a<d;++a) {
                        double ay=actual[(i*q+t)*d+a]/(anorm>0?std::sqrt(anorm):1);
                        double v=ay-sum[a]/norm; error+=v*v;
                    }
                    result[t*3+z]+=error/(anorm>0?1:1e-12)/n;
                    if(++z==3) break;
                }
            }
        }
    }
}
