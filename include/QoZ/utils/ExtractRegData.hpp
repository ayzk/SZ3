
#ifndef SZ_EXTRACTREG_HPP
#define SZ_EXTRACTREG_HPP
namespace QoZ {
	template<class T, uint N>
    inline void extract_lorenzoreg_2d(T *data, std::vector<double> & extracted_A, std::vector<double> & extracted_b, std::vector<size_t> &dims,size_t level,size_t step){
    	assert(dims.size() == N);
    	size_t dimx=dims[0];
    	size_t dimy=dims[1];
    	size_t x_offset=dimy;
    	size_t cur_idx;
    	//size_t y_offset=1;
    	for (size_t i=level;i<dimx;i+=step){
    		for(size_t j=level;j<dimy;j+=step){
    			cur_idx=i*x_offset+j;
    			extracted_b.push_back(data[cur_idx]);
    			for(size_t ii=i-level;ii<=i;ii++){
    				for(size_t jj=j-level;jj<=j;jj++){
    					if(ii==i and jj==j)
    						break;
    				
	    				cur_idx=ii*x_offset+jj;
	    				extracted_A.push_back(data[cur_idx]);
    				}
    			}
    		}
    	}


    }

    template<class T, uint N>
    inline void extract_lorenzoreg_3d(T *data, std::vector<double> & extracted_A, std::vector<double> & extracted_b, std::vector<size_t> &dims,size_t level,size_t step){
    	assert(dims.size() == N);
    	size_t dimx=dims[0];
    	size_t dimy=dims[1];
    	size_t dimz=dims[2];
    	size_t y_offset=dimz;
    	size_t x_offset=y_offset*dimy;
    	size_t cur_idx;
    	//size_t y_offset=1;
    	for (size_t i=level;i<dimx;i+=step){
    		for(size_t j=level;j<dimy;j+=step){
    			for(size_t k=level;k<dimz;k+=step){
	    			cur_idx=i*x_offset+j*y_offset+k;
	    			extracted_b.push_back(data[cur_idx]);
	    			for(size_t ii=i-level;ii<=i;ii++){
	    				for(size_t jj=j-level;jj<=j;jj++){
	    					for(size_t kk=k-level;kk<=k;kk++){
		    					if(ii==i and jj==j and kk==k)
		    						break;
	    				
			    				cur_idx=ii*x_offset+jj*y_offset+kk;
			    				extracted_A.push_back(data[cur_idx]);
			    			}
	    				}
	    			}
    			}
    		}
    	}


    }

}
#endif



