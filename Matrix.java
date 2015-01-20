/*
 * Matrix class to perform Matrix operations
 * By Avi Chad-Friedman 
 * 10/28/14
 */
public class Matrix {
	double[][] data;
	public int rows;
	public int columns;
	
	public Matrix(int rows, int columns){
		this.rows = rows;
		this.columns = columns;
		data = new double[rows][columns];
	}
	public Matrix(Matrix A){
		this.rows = A.rows;
		this.columns = A.columns;
		for(int i=0; i<rows; ++i){
			for(int j=0; j<columns; ++j){
				this.data[i][j] = A.data[i][j];
			}
		}
	}
	
	public Matrix scalarMult(double a){
		Matrix A = this;
		for(int i=0; i<rows; i++){
			for(int j=0; j<columns; j++){
				A.data[i][j] *= a;
			}
		}
		return A;
	}
	public double norm(){
		Matrix A = this;
		if(A.columns != 1)
			throw new RuntimeException("Illegal vector dimensions");
        double norm = 0.0; 
		for(int i = 0; i<A.rows; i++){
        	 norm += A.data[i][0]*A.data[i][0];
         }
		norm = Math.pow(norm, 0.5);
		return norm;
	}
	public Matrix plus(double a){
		Matrix A = this;
		for(int i=0; i<rows; i++){
			for(int j=0; j<columns; j++){
				A.data[i][j] += a;
			}
		}
		return A;
	}
	public Matrix minus(double a){
		Matrix A = this;
		for(int i=0; i<rows; i++){
			for(int j=0; j<columns; j++){
				A.data[i][j] -= a;
			}
		}
		return A;
	}
	public void swapRows(int i, int j){
		Matrix A = this;
		double[] temp = A.data[i];
		A.data[i] = A.data[j];
		A.data[j] = temp;
	}
    public Matrix transpose() {
        Matrix A = new Matrix(columns, rows);
        for (int i = 0; i < rows; i++){
            for (int j = 0; j < columns; j++){
                A.data[j][i] = this.data[i][j];
            }
        }
        return A;
    }
    public Matrix plus(Matrix B){
       Matrix A = this;
       if(A.rows != B.rows || A.columns != B.columns){
    	   throw new RuntimeException("Unequal matrix dimensions");
       }
	   for (int i = 0; i < rows; i++){
          for (int j = 0; j < columns; j++){
        	  A.data[i][j] += B.data[i][j];
          }
	   }
	   return A;	  
    }
    public Matrix minus(Matrix B){
        Matrix A = this;
        if(A.rows != B.rows || A.columns != B.columns){
     	   throw new RuntimeException("Unequal matrix dimensions");
        }
 	   for (int i = 0; i < rows; i++){
           for (int j = 0; j < columns; j++){
         	  A.data[i][j] -= B.data[i][j];
           }
 	   }
 	   return A;	  
     }
    public Matrix matrixMult(Matrix B){
    	Matrix A = this;
        if (A.columns != B.rows) throw new RuntimeException("Illegal matrix dimensions.");
        Matrix C = new Matrix(A.rows, B.columns);
        for (int i = 0; i < C.rows; i++)
            for (int j = 0; j < C.columns; j++)
                for (int k = 0; k < A.columns; k++)
                    C.data[i][j] += (A.data[i][k] * B.data[k][j]);
        return C;
    }
    
    //Returns the inverse
    public Matrix invert(){
    	Matrix A = this;
    	Matrix inverse = new Matrix(A.rows, A.columns);
    	//initialize inverse as identity matrix
    	for(int i = 0; i < inverse.rows; i++){
    		inverse.data[i][i] = 1.0;
    	}
        // Gaussian elimination with partial pivoting
        for (int i = 0; i < rows; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < rows; j++)
                if (Math.abs(A.data[j][i]) > Math.abs(A.data[max][i]))
                    max = j;
            A.swapRows(i, max);
            inverse.swapRows(i, max);
            if (A.data[i][i] == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within A
            for (int j = i + 1; j < rows; j++) {
                double m = A.data[j][i] / A.data[i][i];
                for (int k = 0; k < columns; k++) {
                    A.data[j][k] -= A.data[i][k] * m;
                    inverse.data[j][k] -= inverse.data[i][k] * m;
                }
                A.data[j][i] = 0.0;
            }
        }

        int i = A.rows - 1;
    	int j = A.columns - 1;
        while(i > 0){
        	double piv = A.data[i][j];
    		int k;
    		for(k = i-1; k >= 0; k--){
    			if(A.data[k][j] != 0){
    				break;
    			}
    		}
    		if(A.data[k][j] != 0){
    			double m = A.data[k][j] / piv;
    			A.data[k][j] = 0;
    			for(int t = A.columns-1; t>=0; -- t){
    				inverse.data[k][t] -= m*inverse.data[i][t];
    			}
    		}
    		i--;
        	j--;
    	}
        for(int k = 0; k<A.rows; k++){
        	for(int s = 0; s<A.columns; s++){
        		inverse.data[k][s] = inverse.data[k][s]/A.data[k][k];
        	}
        	A.data[k][k] = 1;
        }
        return inverse;
    }
 // return x = A^-1 b, assuming A is square and has full rank
    public Matrix solve(Matrix rhs) {
        if (rows != columns || rhs.rows != columns || rhs.columns != 1)
            throw new RuntimeException("Illegal matrix dimensions.");

        // create copies of the data
        Matrix A = new Matrix(this);
        Matrix b = new Matrix(rhs);

        // Gaussian elimination with partial pivoting
        for (int i = 0; i < columns; i++) {

            // find pivot row and swap
            int max = i;
            for (int j = i + 1; j < columns; j++)
                if (Math.abs(A.data[j][i]) > Math.abs(A.data[max][i]))
                    max = j;
            A.swapRows(i, max);
            b.swapRows(i, max);

            // singular
            if (A.data[i][i] == 0.0) throw new RuntimeException("Matrix is singular.");

            // pivot within b
            for (int j = i + 1; j < columns; j++)
                b.data[j][0] -= b.data[i][0] * A.data[j][i] / A.data[i][i];

            // pivot within A
            for (int j = i + 1; j < columns; j++) {
                double m = A.data[j][i] / A.data[i][i];
                for (int k = i+1; k < columns; k++) {
                    A.data[j][k] -= A.data[i][k] * m;
                }
                A.data[j][i] = 0.0;
            }
        }
        // back substitution
        Matrix x = new Matrix(columns, 1);
        for (int j = columns - 1; j >= 0; j--) {
            double t = 0.0;
            for (int k = j + 1; k < columns; k++)
                t += A.data[j][k] * x.data[k][0];
            x.data[j][0] = (b.data[j][0] - t) / A.data[j][j];
        }
        return x;
    }
    
    //Perform Gaussian elimination
    public Matrix gaussianEliminate(){
    	Matrix A = this;
    	
    	int i = 0;
    	int j = 0;
    	while(i<A.rows){
    		boolean pivotExists = false;
    		while(j<A.columns && !pivotExists){
    			System.out.println(A.data[i][j]);
    			if(A.data[i][j] != 0){
    				pivotExists = true;
    			}
    			//check for swaps
    			else{
    				int maxRow = i;
    				double maxValue = 0;
    				for(int k = i+1; k<rows; ++k){
    					double curValue = A.data[i][j];
    					if(curValue>maxValue){
    						maxValue = curValue;
    						maxRow = k;
    					}
    				}
    				if(maxRow != i){
    					A.swapRows(maxRow, i);
    					pivotExists = true;
    				}
    				else{
    					j++;
    				}
    		}
    			
    		}
    		if(pivotExists){
    			for(int t = i + 1; t<A.rows; ++t){
    				for(int s=j+1; s<A.columns; ++s){

    					A.data[t][s] = A.data[t][s] - A.data[i][s]*(A.data[t][j]/A.data[i][j]);
    				}
    				A.data[t][j] = 0;
    			}
    		}
    		i++;
    		j++;
    		
    	}
    	return A;
    }
    
    //Matrix must already have been gaussian eliminated
    public Matrix reducedRowForm(){
    	Matrix A = this;
    	int i = A.rows - 1;
    	int j = A.columns - 2;
    	while(i >= 0){
    		int k = j-1;
    		while(k>=0){
    			if(A.data[i][k] != 0){
    				j = k;
    			}
    			k--;
    		}
    		//Make elements above the pivots 0
    		if(A.data[i][j] != 0){
    			for(int t = i-1; t>=0; --t){
    				for(int s = 0; s<A.columns; ++s){
    					if(s != j){
    						A.data[t][s] = A.data[t][s]-A.data[i][s]*(A.data[t][j]/A.data[i][j]);
    					}
    				}
    				A.data[t][j] = 0;
    			}
    		
    		for(k=j+1; k<A.columns; ++k){
    			A.data[i][k] = A.data[i][k]/A.data[i][j];
    		}
    		A.data[i][j] = 1;
    		}
    		i--;
        	j--;
    	}
    	return A;
    }
    
    
    public void printAnswers(){
    	Matrix A = this;
    	boolean solutionExists = true;
    	boolean done = false;
    	int i = 0;
    	while(!done && i<A.rows){
    		boolean allZeros = true;
    		for(int j=0; j<A.columns-1; ++j){
    			if(A.data[i][j] != 0){
    				allZeros = false;
    			}
    		}
    		//If the row of R is all zeros but d isn't
    		if(allZeros && A.data[i][A.columns-1] != 0){
    			solutionExists = false;
    			System.out.println("No solution.");
    		}
    		else if(allZeros && A.data[i][A.columns-1] == 0){
    			done = true;
    			System.out.println("Infinite Solutions");
    		}
    		else if(A.rows<(A.columns-1)){
    			done = true;
    			System.out.println("Infinite solutions");
    		}
    		i++;
    	}
    	if(!done){
    		System.out.println("Unique Solution");
    	}
    	if(solutionExists){
    		Matrix particular = new Matrix(A.columns-1, 1);
    		boolean[] freeMarker = new boolean[A.columns-1];
			int n = 0;
    		for(i=0; i<A.rows; ++i){
    			boolean pivotFound = false;
        		//int free2 = 0;
    			for(int j=0; j<A.columns-1; ++j){
    				if(A.data[i][j] != 0){
    					//If the pivot is found, add to particular
    					if(!pivotFound){
    						pivotFound = true;
        					particular.data[j][0] = A.data[i][A.columns-1];
        					freeMarker[j] = false;
    					}
    					//Only want to create one set of special solutions
    					else{
    						freeMarker[j] = true;
    					}
    				
    				}
    				
    			}
    		}
    		System.out.println("Particular Solution");
    		
    		particular.printMatrix();
    		for(int j = 0; j<A.columns-1; j++){
    			//Create a special solution for each free variable
    			if(freeMarker[j]){
    				Matrix special = new Matrix(A.columns-1, 1);
    				special.data[j][0] = 1;
    				for(int k = 0; k<A.columns-1; k++){
    					//Set the non-free variables of the special solutions to the negated values of column
    					if(!freeMarker[k]){
    						for(int z = 0; z<A.rows; z++){
    							if(A.data[k][z] == 1){
    								special.data[k][0] = -1*A.data[z][j];
    							}
    						}
    					}
    				}
    				System.out.println("Special Solution");
            		special.printMatrix();
    			}
    			
    		}
    		
    	}
    }
    
  public void printMatrix(){
	  Matrix A = this;
	  for(int i = 0; i<A.rows; i++){
		  for(int j = 0; j<A.columns; j++){
			  System.out.print(A.data[i][j]+" ");
		  }
		  System.out.println();
	  }
  }
    
}
