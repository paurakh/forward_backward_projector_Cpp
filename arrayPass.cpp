	#include<iostream>
using namespace std;
void printMatrix(int matrix[4][4]){
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			cout<<matrix[i][j]<<" ";   
		}   
		cout<<endl;
	}

}

void printarray( int **array ){
    cout << array[0][0];
}

int main(){


    // initialize 2d pointer to pointer (~array)
	int row = 10, col = 5;
	int **board;
	board = new int*[row]; // dynamic array (size 10) of pointers to int
	for (int i = 0; i < row; ++i) {
		board[i] = new int[col];
	  // each i-th pointer is now pointing to dynamic array (size 10) of actual int values
	}


    // example of assigning variables
	for (int i = 0; i < row; ++i) {   // for each row
	  for (int j = 0; j < col; ++j) { // for each column
	  	board[i][j] = i*j;
	  	cout << board[i][j] << " ";   // sample print each element followed by a space

	  }
	  cout << endl; // print line break
	}

	printarray(board);

	// deallocate the memory
	for (int i = 0; i < col; i++) {
		delete[] board[i];
	}
	delete[] board;
	return 0;
}