/**
 * Created by chenyi(wugu) on 2017/3/3.
 */
//line行 row列
public class LinearAlgebra {
    //判断能否进行矩阵乘法
    private static boolean canMultiply(Matrix left,Matrix right) {
        return left.row == right.line;
    }
    //矩阵乘法
     private static Matrix multiply(Matrix left,Matrix right) {
        if (canMultiply(left,right)) {
            Matrix result = new Matrix(new double[left.line][right.row]);
            for (int l = 0; l < left.line; l++) {
                for (int r = 0; r < right.row; r++) {
                    double sum = 0;
                    for (int u = 0; u < left.row; u++) {
                        sum = sum + left.matrix[l][u] * right.matrix[u][r];
                    }
                    result.matrix[l][r] = sum;
                }
            }
            return result;
        } else {
            return null;
        }
    }
    //判断能否进行矩阵加法
    private static boolean canAdd(Matrix left,Matrix right) {
        return (left.line == right.line && left.row == right.row);
    }
    //矩阵加法
    private static Matrix add(Matrix left,Matrix right) {
        if (canAdd(left,right)) {
            Matrix result = new Matrix(new double[left.line][left.row]);
            for (int l = 0; l < left.line; l++) {
                for (int r = 0; r < left.row; r++) {
                    result.matrix[l][r] = left.matrix[l][r] + right.matrix[l][r];
                }
            }
            return result;
        } else {
            return null;
        }
    }
    //矩阵转置
    private static Matrix transpose(Matrix matrix) {
        Matrix result = new Matrix(new double[matrix.row][matrix.line]);
        for (int l = 0; l < matrix.line; l++) {
            for (int r = 0; r < matrix.row; r++) {
                result.matrix[r][l] = matrix.matrix[l][r];
            }
        }
        return result;
    }
    //矩阵减法
    private static Matrix subtract(Matrix left,Matrix right) {
        if (canAdd(left,right)) {
            Matrix result = new Matrix(new double[left.line][left.row]);
            for (int l = 0; l < left.line; l++) {
                for (int r = 0; r < left.row; r++) {
                    result.matrix[l][r] = left.matrix[l][r] - right.matrix[l][r];
                }
            }
            return result;
        } else {
            return null;
        }
    }
    //判断是否为方形矩阵 (能否求行列式)
    private static boolean isSquareMatrix(Matrix matrix) {
        return matrix.line == matrix.row;
    }
    //求行列式
    private static double getDetaminate(Matrix matrix) {
        int i, j, m, n, s, t, k = 1;
        int N = matrix.line;
        double f = 1, c, x, det = 0;
        for (i = 0, j = 0; i < N && j < N; i++, j++) {
            if (matrix.matrix[i][j] == 0) {
                for (m = i; m < matrix.line && matrix.matrix[m][j] == 0; m++) ;
                if (m == N) {
                    det = 0;
                    return det;
                } else
                    for (n = j; n < N; n++) {
                        c = matrix.matrix[i][n];
                        matrix.matrix[i][n] = matrix.matrix[m][n];
                        matrix.matrix[m][n] = c;
                    }
                k *= (-1);
            }
            for (s = N - 1; s > i; s--) {
                x = matrix.matrix[s][j];
                for (t = j; t < N; t++)
                    matrix.matrix[s][t] -= matrix.matrix[i][t] * (x / matrix.matrix[i][j]);
            }
        }
        for (i = 0; i < N; i++)
            f *= matrix.matrix[i][i];
        det = k * f;
        return det;
    }

    //求伴随矩阵
    private static Matrix getAdjoint(Matrix matrix) {
        double[][][] matrix_z = new double[matrix.line * matrix.row][matrix.line - 1][matrix.row - 1];

        for (int k = 0, z = 0; k < matrix.line && z < matrix.line * matrix.row; k++) {
            for (int K = 0; K <matrix.row; K++) {
                for (int i5 = 0, m5 = 0; i5 < matrix.line && m5 < matrix.line - 1; i5++) {
                    for (int p = 0, M = 0; p < matrix.row && M < matrix.row - 1; p++) {
                        if (p != K) {
                            if ((i5 + p) % 2 == 0) {
                                matrix_z[z][M][m5] = matrix.matrix[i5][p];
                                M++;
                            } else {
                                matrix_z[z][M][m5] = -matrix.matrix[i5][p];
                                M++;
                            }
                        }
                    }
                    if (i5 != k)
                        m5++;
                }
                z++;
            }
        }
        double[][] resultArr = new double[matrix.line][matrix.row];
        for (int count = 0, l = 0; l < matrix.line; l++) {
            for (int r = 0; r < matrix.row; r++, count++) {
                Matrix matrixB = new Matrix(matrix_z[count]);
                resultArr[l][r] = getDetaminate(matrixB);
            }
        }
        Matrix result = new Matrix(resultArr);
        return result;
    }
    //求逆矩阵
    private static Matrix getInverseMatrix(Matrix matrix) {
        Matrix result = getAdjoint(matrix);
        double det = getDetaminate(matrix);
        for (int l = 0; l < matrix.line; l++) {
            for (int r = 0; r < matrix.row; r++) {
                result.matrix[l][r] = result.matrix[l][r] / det;
            }
        }
        return result;
    }
    //矩阵正交化
    private static Matrix getOrthogonalization(Matrix matrix) {
        Matrix result = new Matrix(new double[matrix.line][matrix.row]);
        for (int i1 = 0; i1 < matrix.line; i1++)//行
        {
            for (int i4 = 0; i4 < matrix.row; i4++)//列
            {
                double sum3 = 0;
                for (int i3 = 0; i3 < i1; i3++) {
                    double sum1 = 0, sum2 = 0;

                    for (int i2 = 0; i2 < matrix.row; i2++) {
                        sum1 += matrix.matrix[i1][i2] * result.matrix[i3][i2];
                        sum2 += result.matrix[i3][i2] * result.matrix[i3][i2];
                    }
                    sum3 += (sum1 / sum2) * result.matrix[i3][i4];
                }
                result.matrix[i1][i4] = matrix.matrix[i1][i4] - sum3;
            }
        }
        return result;
    }
    //矩阵正交单位化
    private static Matrix getOrthogonalUnit(Matrix matrix) {
        Matrix result = getOrthogonalization(matrix);
        return getUnit(matrix);
    }
    private static Matrix getUnit(Matrix matrix){
        Matrix result = new Matrix(new double[matrix.line][matrix.row]);
        for(int l = 0;l<matrix.line;l++){
            double sum = 0;
            for(int r = 0;r<matrix.row;r++){
                sum += matrix.matrix[l][r]* matrix.matrix[l][r];
            }
            for(int r = 0;r<matrix.row;r++){
                result.matrix[l][r] = matrix.matrix[l][r]/Math.sqrt(sum);
            }
        }
        return result;
    }


    //test
    public static void main(String[] args) {
        double[][] c = {{1,2,4},{2,-2,2},{4,2,1}};
        Matrix a = new Matrix(c);
        Matrix b = new Matrix(c);
        for(int i =1;i<10;i++){
            b = multiply(b,a);
            System.out.println(i+":");
            b.print();
        }
        b.print();
    }
}

