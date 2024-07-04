#include "algebra.h"
#include <stdio.h>
#include <math.h>

Matrix create_matrix(int row, int col)
{
    Matrix m;
    m.rows = row;
    m.cols = col;
    return m;
}

Matrix add_matrix(Matrix a, Matrix b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0); 
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] + b.data[i][j];
        }
    }
    return result;
}
Matrix sub_matrix(Matrix a, Matrix b) {
    if (a.rows != b.rows || a.cols != b.cols) {
        printf("Error: Matrix a and b must have the same rows and cols.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] - b.data[i][j];
        }
    }
    return result;
}
Matrix mul_matrix(Matrix a, Matrix b) {
    if (a.cols != b.rows) {
        printf("Error: The number of cols of matrix a must be equal to the number of rows of matrix b.\n");
        return create_matrix(0, 0);
    }
    Matrix result = create_matrix(a.rows, b.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < b.cols; j++) {
            result.data[i][j] = 0;
            for (int k = 0; k < a.cols; k++) {
                result.data[i][j] += a.data[i][k] * b.data[k][j];
            }
        }
    }
    return result;
}
Matrix scale_matrix(Matrix a, double k) {
    Matrix result = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[i][j] = a.data[i][j] * k;
        }
    }
    return result;
}

Matrix transpose_matrix(Matrix a) {
    Matrix result = create_matrix(a.cols, a.rows);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            result.data[j][i] = a.data[i][j];
        }
    }
    return result;
}

double det_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }

    if (a.rows == 1)return a.data[0][0];
    
    double result = 0.0;

    for (int i = 0; i < a.rows; i++) {
        Matrix submatrix = create_matrix(a.rows - 1, a.rows - 1);
        for (int row = 1; row < a.rows; row++) {
            int colCount = 0;
            for (int col = 0; col < a.rows; col++) {
                if (col == i) {
                    continue;
                }
                submatrix.data[row - 1][colCount] = a.data[row][col];
                colCount++;
            }
        }
        double sub_det = det_matrix(submatrix);
        result += (i % 2 == 0 ? 1 : -1) * a.data[0][i] * sub_det;
    }
    return result;
}

Matrix inv_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return create_matrix(0, 0);
    }
    double det = det_matrix(a);
    if (det == 0) {
        printf("Error: The matrix is singular.\n");
        return create_matrix(0, 0);
    }
    Matrix adj = create_matrix(a.rows, a.cols);
    for (int i = 0; i < a.rows; i++) {
        for (int j = 0; j < a.cols; j++) {
            Matrix sub = create_matrix(a.rows - 1, a.cols - 1);
            int sub_i = 0, sub_j = 0;
            for (int row = 0; row < a.rows; row++) {
                if (row == i) continue;
                sub_j = 0;
                for (int col = 0; col < a.cols; col++) {
                    if (col == j) continue;
                    sub.data[sub_i][sub_j] = a.data[row][col];
                    sub_j++;
                }
                sub_i++;
            }
            double sub_det = det_matrix(sub);
            adj.data[j][i] = ((i + j) % 2 == 0 ? 1 : -1) * sub_det;
            if (adj.data[j][i] == -0) adj.data[j][i] = 0;
        }
    }
    Matrix inv = scale_matrix(adj, 1.0 / det);
    return inv;
}

int rank_matrix(Matrix a) {
    int m = a.rows;
    int n = a.cols;
    int rank = (m < n) ? m : n;

    for (int row = 0; row < rank; row++) {
        if (a.data[row][row] == 0) {
            int reduce = 1;
            for (int k = row + 1; k < m; k++) {
                if (a.data[k][row] != 0) {
                    for (int j = 0; j < n; j++) {
                        double temp = a.data[row][j];
                        a.data[row][j] = a.data[k][j];
                        a.data[k][j] = temp;
                    }
                    reduce = 0;
                    break;
                }
            }
            if (reduce) break;}
        for (int r = row + 1; r < m; r++) {
            if (a.data[r][row] != 0) {
                double mult = a.data[r][row] / a.data[row][row];
                for (int col = 0; col < n; col++) {
                    a.data[r][col] -= mult * a.data[row][col];
                }
            }
        }
    }

        rank = m;
    for (int i = m - 1; i >= 0; i--) {
        int index = 0;
        for (int j = 0; j < n; j++){
            if(fabs(a.data[i][j])>0.000000001){
                index = 1;
                break;
            }
        }
        if (index == 0) rank--;
    }
    return rank;
}

double trace_matrix(Matrix a) {
    if (a.rows != a.cols) {
        printf("Error: The matrix must be a square matrix.\n");
        return 0;
    }
    double trace = 0;
    for (int i = 0; i < a.rows; i++) {
        trace += a.data[i][i];
    }
    return trace;
}

void print_matrix(Matrix a)
{
    for (int i = 0; i < a.rows; i++)
    {
        for (int j = 0; j < a.cols; j++)
        {
            // 按行打印，每个元素占8个字符的宽度，小数点后保留2位，左对齐
            printf("%-8.2f", a.data[i][j]);
        }
        printf("\n");
    }
}