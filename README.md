# SMR

## Requirement  

* Eigen
* libomp or libgomp (OpenMP is implemented in GCC, you need to modify your gcc except that you want to link libomp staticly.)    
* libz  

## install

### install requirments  

put Eigen into your $C_INCLUDE_PATH (such as /usr/include, /usr/local/include)  
check if libomp(if there just libgomp, please make a link) and libz was installed in you operation system. some Linux distribution maybe need install these lib from repository manually.  
### compile

In the top level of source   
`make`  

## Usage  

smr --help  
smr --version  
smr smr_main [\<args\>]  
