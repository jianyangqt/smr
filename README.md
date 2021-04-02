# SMR

## Requirement  

* Eigen
* libomp or libgomp (OpenMP is implemented in GCC)  
* libz  

## install

### install requirments  
put Eigen into your $C_INCLUDE_PATH (such as /usr/include, /usr/local/include)  
check if libomp(if there just libgomp, please make a link) and libz was installed in you operation system. some Linux distribution maybe need install these lib from repository manually.  

### compile

In the top level of source   
`./configure`  
`make`  
`make install`  
