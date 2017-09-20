/**
 * Created by Lenovo on 2017/3/5.
 */
class Matrix {

    double[][] matrix;
    int line;
    int row;
     Matrix(double[][] matrix) {
        this.matrix = matrix;
        this.line = matrix.length;
        this.row = matrix[0].length;
    }
    //打印矩阵
    public void print() {
        for (int l = 0; l < line; l++) {
            for (int r = 0; r < row; r++) {
                System.out.format("%.3f   ",matrix[l][r]);
            }
            System.out.println();
        }
    }
}
