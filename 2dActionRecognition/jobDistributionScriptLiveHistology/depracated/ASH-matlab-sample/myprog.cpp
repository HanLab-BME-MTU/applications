#include <iostream>
#include <cstdlib>
using namespace std;

void help(){
	cout << "myprog m n" << endl;
}

void foo(int m, int n){
	const int N = 100000000;//0.1G
	int result = 0;
	for(int i = 0; i < N; i++){
		for(int j = 0; j < 300; j++){
			result += m * n;
		}
	}
}

int main(int argc, char** argv)
{
	if(argc != 3){
		help();
		return 0;
	}
	
	int m = atoi(argv[1]);
	int n = atoi(argv[2]);
	
	foo(m, n);
	
	cout << "myprog finished." << endl;
	return 0;
}
