#include<iostream>
#include<fstream>
#include<math.h>
#include<iomanip>
#include<cstdlib>
using namespace std;
double rounding(double num, int index)
{
    bool isNegative = false; // whether is negative number or not
	
    if(num < 0) // if this number is negative, then convert to positive number
    {
        isNegative = true;	
        num = -num;
    }
	
    if(index >= 0)
    {
        int multiplier = 1;
        for (int h = 1; h <= index;h++)
            multiplier *= 10;
        num = (int)(num * multiplier + 0.5) / (multiplier * 1.0);
    }
	
    if(isNegative) // if this number is negative, then convert to negative number
    {
        num = -num;
    }
    return num;
}
double calculate_q(int column_index,double A[][11],double q[][11]){  //轉成單位向量
    float c = 0;
    for (int i = 0; i < 11;i++)
    {
        c += pow(A[i][column_index],2);
    }
    c = pow(c, 0.5);
    if(c<0.0001)
        return 0;
    for (int i = 0; i < 11;i++)
        {   
            q[i][column_index] = A[i][column_index] / c;
        }
    return 0;
}
void calculate_A(int column_index, double a[][11], double q[][11],double A[][11]){
    double array[11] = {0};
            int m = 1;
    for (int i = 1; i < column_index;i++)
        {
            double c = 0; //記錄內積結果
            for (int j = 1; j < 11;j++)
            {
                c += a[j][column_index] * q[j][m];
            }
            for (int j = 1; j < 11;j++)
            {
                array[j] += c *q[j][m];
            }
                m++;
        }
        for (int i = 1; i < 11; i++)
        {
            A[i][column_index] = a[i][column_index] - array[i];
        }
    calculate_q(column_index, A, q);//算出q向量(把A向量轉為單位向量)
}
void calculate_r(double a[][11], double q[][11],double r[][11]){
    for (int i = 1; i < 11;i++)
        for (int j = i; j < 11;j++)
        {
            double c = 0;
            for (int k = 1; k < 11;k++)
                c += a[k][j] * q[k][i];
            r[i][j] = c;
        }
}
void multiply(double one[][11],double two[][11],double ans[][11])//矩陣乘法
{
    for (int i = 0; i < 11;i++)
        for (int j = 0; j < 11;j++)
            ans[i][j] = 0;
    for (int i = 1; i < 11; i++)
        for (int j = 1; j < 11; j++)
            for (int k = 1; k < 11; k++)
                ans[i][j] += one[i][k] * two[k][j];
}
void elementary_2(int m,int n,double (&a)[11][11],double (&u)[11][21],double V[][11])
{
    double c ;
    int v = 1;
    bool result = false;
    int p[18][11][11];
    double a_row_exchanged[11][11] = {{0}};
    int not_zero;
    for (int y = 1; y <= 10;y++)
        for (int q = 1; q <= 20;q++)
            u[y][q] = 0;
    for (int y = 1; y <= 10;y++)
        for (int q = 1; q <= 10;q++)
            u[y][q] = a[y][q]; //把a的內容複製到u的左側
    for (int y = 1; y <= 10;y++)
            u[y][y+10] = 1;    //把u的右側初始化為單位矩陣
    for (int i = 1; i < m&&v<=n;i++) //每個row都執行一次往下消去
    {
        result = false;
        for (int f = i+1; f <= m;f++)
                if(u[f][v]!=0)
                    result = true;
        if (u[i][v] == 0&& (result==true))//pivot在下面某行
            {
                for (int s = i+1; s <= m;s++)
                    if(u[s][i]!=0)
                        not_zero = s;
                for (int l = 1; l <= 20; l++)
                    {
                        int temp;
                        temp = u[i][l];
                        u[i][l] = u[not_zero][l];
                        u[not_zero][l] = temp;
                    }//交換整個row
                    p[i][i][i] = 0;
                    p[i][not_zero][not_zero] = 0;
                    p[i][i][not_zero] = 1;
                    p[i][not_zero][i] = 1;
            }//把i row與下面某一列起始不為零的row交換
        for (int k = i + 1; k <= m; k++)//對第k行消去
        {
            if(result==false&&u[i][v] == 0)
                {
                    v++;
                    k--;
                    continue;
                }
            if(u[k][i]==0)
                continue;
            c = -(u[k][v] / u[i][v]);
            for (int y = v; y <= 20; y++) //處理第k行
                u[k][y] += u[i][y] * c;
        }
        v++;
    }
        for (int i = 1; i <9 ;i++)
            for (int j = 1; j <= 10;j++)
                for (int k = 1; k <= 10;k++)
                        {
                        for (int q = 1; q <= 10;q++)
                            p[9+i][j][k] += (p[8+i][j][q] * p[9-i][q][k]);
                        }
        for (int i = 1; i <= 10;i++)
                for (int j = 1; j <= 10;j++)
                    {
                        for (int q = 1; q <= 10;q++)
                            a_row_exchanged[i][j] += p[17][i][q] * a[q][j];
                    }
        for (int y = 1; y <= 10;y++)
            for (int q = 1; q <= 10;q++)
                u[y][q] = a[y][q]; //把a的內容複製到u的左側
        for (int y = 1; y <= 10;y++)
            u[y][y+10] = 1;    //把u的右側初始化為單位矩陣
        v = 0;
        for (int i = 1; i < m&&v<=n;i++) //每個row都執行一次往下消去
    {
        for (int k = i + 1; k <= m; k++)//對第k行消去
        {
            if(u[k][v]==0)
                continue;
            else
                c = -(u[k][v] / u[i][v]);
            result = false;
            for (int f = i+1; f <= m;f++)
                if(u[f][v]!=0)
                    result = true;
            if(result==false&&u[i][v] == 0)
                {
                    v++;
                    k--;
                    continue;
                }
            if (u[i][v] == 0&& (result==true))
                return ; //如果某列的第i個element已經被消成0，代表需要p矩陣
            for (int p = v; p <= 20; p++) //處理i的下面一行
                u[k][p] += u[i][p] * c;
        }
        v++;
    }
    return;
}
void elementary_1(int m,int n,double (&a)[11][11],double (&u)[11][21],double V[][11])
{
    double c;
    int v = 1;
    bool result = false;
    for (int y = 1; y <= 10;y++)
        for (int q = 1; q <= 20;q++)
            u[y][q] = 0;
    for (int y = 1; y <= 10;y++)
        for (int q = 1; q <= 10;q++)
            u[y][q] = a[y][q]; //把a的內容複製到u的左側
    for (int y = 1; y <= 10;y++)
            u[y][y+10] = 1;    //把u的右側初始化為單位矩陣
    for (int i = 1; i <m&&v<=n;i++) //每個row都執行一次往下消去
    {
        for (int k = i + 1; k <= n; k++)//對第k列消去
        {
            if(u[k][v]==0)
                continue;
            else
                c = -(u[k][v] / u[i][v]);
            result = false;
            for (int f = i+1; f <= m;f++)
                if(u[f][v]!=0)
                    result = true;
            if(result==false&&u[i][v] == 0)//整行都是0
                {
                    v++;
                    k--;//從原本那行繼續往下消
                    continue;
                }
            if(u[i][v] == 0&& (result==true))//pivot在別行
                {
                    elementary_2(m, n, a, u,V);
                    return;
                } //如果某列的第i個element已經被消成0，代表需要p矩陣
            for (int p = v; p <= 20; p++) //處理第k列
                u[k][p] += u[i][p] * c;
        }
        v++;
    }
    return;
}
void QR_algorithm(int m,int n,int r,double ATA[][11],double QV[][11],double S[][11])
{
    int width_length;
    if(m==r)
        width_length = m;
    else
        width_length = n;
    int column_of_eigenvector = 1;
    float lambda[11] = {0}; //用來存eigen value
    double AA[11][11] = {{0}};
    for (int i = 0; i < 11;i++)
            for (int j = 0; j < 11;j++)
                AA[i][j] = ATA[i][j];
    double ATAT[11][11] = {{0}}, Q[11][11] = {{0}}, R[11][11] = {{0}}, A[11][11] = {{0}};
    for (int k = 0; k < 60;k++)
    {
        for (int i = 0; i < 11;i++)
            for (int j = 0; j < 11;j++)
                Q[i][j] = 0;
        for (int i = 1; i < 11;i++) calculate_A(i, AA, Q, A);calculate_r(AA,Q,R);//QR分解
        multiply(R, Q, AA);       //將新的值存進AA
    }
    for (int i = 1,k=1; i <= r;i++)
        {
            for (int j = 1; j <= r;j++)
                if(AA[i][i] == lambda[j])
                        break;
            lambda[k] = AA[i][i];
            k++;
        }
        for (int i = 1; i <= r;i++)
            S[i][i] = pow(lambda[i], 0.5);//完成SIGMA矩陣
        double null_space[11][21] = {{0}}; //用於找nullspace
    for (int k = 1; k <= width_length;k++){
        for (int i = 0; i < 11; i++)
                for (int j = 0; j < 11; j++)
                    ATAT[j][i] = ATA[i][j]; //建立轉置矩陣
        for (int i = 1; i <= width_length; i++)
            ATAT[i][i] -= lambda[k];
        elementary_1(width_length, width_length, ATAT, null_space, QV);
        for (int i = 0; i < 11;i++)
            for (int j = 0; j < 21;j++)
                if((null_space[i][j]<0.01&&null_space[i][j]>0)||(null_space[i][j]>-0.01&&null_space[i][j]<0))
                    null_space[i][j] = 0;
        for (int i = 1; i <=width_length;i++)//檢查是否整列為0
            {
                bool is_eigenvector = true;
                for (int j = 1; j < 11;j++)
                    if(null_space[i][j]!=0)
                        is_eigenvector = false;
                if(is_eigenvector)
                    {
                        for (int m = 1; m <= width_length;m++)
                            QV[m][column_of_eigenvector] = null_space[i][m + 10];
                        column_of_eigenvector++;
                    }
            }
                if (lambda[k] == 0)
                    break;
        }
}  
int main(void)
{
    int case_number;
    ifstream inFile("input.txt", ios::in);
    if(!inFile){
        cerr << "Failed opening" << endl;
        exit(1);
    }
    ofstream outFile("output.txt", ios::out);
    if(!outFile){
        cerr << "Failed opening" << endl;
        exit(1);
    }
    inFile >> case_number;
    outFile << case_number << endl;
    while(case_number>0)
    {
        double U[11][11] = {{0}};
        double V[11][11] = {{0}};
        double QV[11][11] = {{0}};
        double QU[11][11] = {{0}};
        double temp[11][11] = {{0}};
        double S[11][11] = {{0}};
        double A[11][11] = {{0}};
        double AT[11][11] = {{0}};
        double ATA[11][11] = {{0}};
        double AAT[11][11] = {{0}};
        int m, n;
        inFile >> m >> n;
        for (int i = 1; i <= m;i++)
            for (int j = 1; j <= n;j++)
                inFile >> A[i][j];//輸入待分解矩陣
        for (int i = 0; i < 11;i++)
            for (int j = 0; j < 11;j++)
                AT[j][i] = A[i][j];//建立轉置矩陣
        multiply(AT, A, ATA);
        QR_algorithm(m,n,n,ATA,V,S);//QR_algorithm輸入(m,n,r,ATA/AAT,V/U,S)
        for (int i = 1; i < 11;i++)
            calculate_A(i, V, QV, temp);
        multiply(A, AT, AAT);
        QR_algorithm(m, n, m, AAT, U, S);
        for (int i = 1; i < 11;i++)
            calculate_A(i, U, QU, temp);
        outFile << m << " " << m << endl;
        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= m;j++)
                {
                    if(rounding(QU[i][j], 2)==0)
                        QU[i][j] *= -1;
                outFile << setiosflags(ios::fixed) << setprecision(2) << rounding(QU[i][j], 2) << " ";
                }
                outFile << endl;
        }
        outFile << m << " " << n << endl;
        for (int i = 1; i <= m; i++)
        {
            for (int j = 1; j <= n;j++)
                {
                    if(rounding(S[i][j], 2)==0)
                        S[i][j] *= -1;
                outFile << setiosflags(ios::fixed) << setprecision(2) << rounding(S[i][j], 2) << " ";
                }
                outFile << endl;
        }
        outFile << n << " " << n << endl;
        for (int i = 1; i <= n; i++)
        {
            for (int j = 1; j <= n;j++)
                {
                    if(rounding(QV[i][j], 2)==0)
                        QV[i][j] *= -1;
                outFile << setiosflags(ios::fixed) << setprecision(2) << rounding(QV[i][j], 2) << " ";
                }
                outFile << endl;
        }
        case_number--;
    }
}