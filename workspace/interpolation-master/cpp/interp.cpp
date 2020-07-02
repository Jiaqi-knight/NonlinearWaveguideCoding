#include "interp.h"

//=====================  BODY  ==========================

//-----------------
// Constructors 
//-----------------

// Constructor (default)
template<class T>
interp<T>::interp(){}

// Constructor (1D)
template<class T>
interp<T>::interp(T xinp[], T winp[], T xvalinp[], int Ninp, int Pinp, int Kinp){
    D  = 1;
    N = Ninp;
    M = 1;
    L = 1;
    P = Pinp;
    Q = 1;
    R = 1;
    K = Kinp;
    x.resize(N);
    xval.resize(P);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[n+k*N]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int p=0; p<P; p++){
        f[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (1D measure performances)
template<class T>
interp<T>::interp(T xinp[], T winp[], T xvalinp[], int Ninp, int Pinp, int Kinp, bool meas){
    D  = 1;
    measure = meas;
    N = Ninp;
    M = 1;
    L = 1;
    P = Pinp;
    Q = 1;
    R = 1;
    K = Kinp;
    x.resize(N);
    xval.resize(P);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[n+k*N]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int p=0; p<P; p++){
        f[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (1D)
template<class T>
interp<T>::interp(vector<T> xinp, vector<vector<T> > winp, vector<T> xvalinp){
    D  = 1;
    N = xinp.size();
    M = 1;
    L = 1;
    P = xvalinp.size();
    Q = 1;
    R = 1;
    K = winp.size();
    x.resize(N);
    xval.resize(P);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[k][n]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int p=0; p<P; p++){
        f[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (1D measure performances)
template<class T>
interp<T>::interp(vector<T> xinp, vector<vector<T> > winp, vector<T> xvalinp, bool meas){
    D  = 1;
    measure = meas;
    N = xinp.size();
    M = 1;
    L = 1;
    P = xvalinp.size();
    Q = 1;
    R = 1;
    K = winp.size();
    x.resize(N);
    xval.resize(P);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[k][n]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int p=0; p<P; p++){
        f[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (2D)
template<class T>
interp<T>::interp(T xinp[], T yinp[], T winp[], T xvalinp[], T yvalinp[], int Ninp, int Minp, int Pinp, int Qinp, int Kinp){
    D  = 2;
    N = Ninp;
    M = Minp;
    L = 1;
    P = Pinp;
    Q = Qinp;
    R = 1;
    K = Kinp;
    x.resize(N);
    y.resize(M);
    xval.resize(P);
    yval.resize(Q);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int m=0; m<M; m++){
        y[m] = yinp[m];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[m+n*M+k*M*N]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int q=0; q<Q; q++){
        yval[q] = yvalinp[q];
    }
    for(unsigned int p=0; p<P; p++){
        f[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (2D measure performances)
template<class T>
interp<T>::interp(T xinp[], T yinp[], T winp[], T xvalinp[], T yvalinp[], int Ninp, int Minp, int Pinp, int Qinp, int Kinp, bool meas){
    D  = 2;
    measure = meas;
    N = Ninp;
    M = Minp;
    L = 1;
    P = Pinp;
    Q = Qinp;
    R = 1;
    K = Kinp;
    x.resize(N);
    y.resize(M);
    xval.resize(P);
    yval.resize(Q);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int m=0; m<M; m++){
        y[m] = yinp[m];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[m+n*M+k*M*N]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int q=0; q<Q; q++){
        yval[q] = yvalinp[q];
    }
    for(unsigned int p=0; p<P; p++){
        f[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (2D)
template<class T>
interp<T>::interp(vector<T> xinp, vector<T> yinp, vector<vector<vector<T> > > winp, vector<T> xvalinp, vector<T> yvalinp){
    D  = 2;
    N = xinp.size();
    M = yinp.size();
    L = 1;
    P = xvalinp.size();
    Q = yvalinp.size();
    R = 1;
    K = winp.size();
    x.resize(N);
    y.resize(M);
    xval.resize(P);
    yval.resize(Q);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int m=0; m<M; m++){
        y[m] = yinp[m];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[k][n][m]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int q=0; q<Q; q++){
        yval[q] = yvalinp[q];
    }
    for(unsigned int p=0; p<P; p++){
        f[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (2D measure performances)
template<class T>
interp<T>::interp(vector<T> xinp, vector<T> yinp, vector<vector<vector<T> > >  winp, vector<T> xvalinp, vector<T> yvalinp, bool meas){
    D  = 2;
    measure = meas;
    N = xinp.size();
    M = yinp.size();
    L = 1;
    P = xvalinp.size();
    Q = yvalinp.size();
    R = 1;
    K = winp.size();
    x.resize(N);
    y.resize(M);
    xval.resize(P);
    yval.resize(Q);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int m=0; m<M; m++){
        y[m] = yinp[m];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[k][n][m]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int q=0; q<Q; q++){
        yval[q] = yvalinp[q];
    }
    for(unsigned int p=0; p<P; p++){
        f[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (3D)
template<class T>
interp<T>::interp(T xinp[], T yinp[], T zinp[], T winp[], T xvalinp[], T yvalinp[], T zvalinp[], int Ninp, int Minp, int Linp, int Pinp, int Qinp, int Rinp){
    D  = 3;
    N = Ninp;
    M = Minp;
    L = Linp;
    P = Pinp;
    Q = Qinp;
    R = Rinp;
    K = Kinp;
    x.resize(N);
    y.resize(M);
    z.resize(L);
    xval.resize(P);
    yval.resize(Q);
    zval.resize(R);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int m=0; m<M; m++){
        y[m] = yinp[m];
    }
    for(unsigned int l=0; l<L; l++){
        z[l] = zinp[l];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[k][n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[k][n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[k][n][m][l] = winp[l+m*L+n*L*M+k*L*M*N]; 
                }
            }
        }
    }
    for(unsigned int n=0; n<N; n++){
        w[n].resize(M);
        for(unsigned int m=0; m<M; m++){
            w[n][m].resize(L);
            for(unsigned int l=0; l<L; l++){
                w[n][m][l] = winp[l+m*L+n*M*L]; 
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int q=0; q<Q; q++){
        yval[q] = yvalinp[q];
    }
    for(unsigned int r=0; r<R; r++){
        zval[r] = zvalinp[r];
    }
    for(unsigned int p=0; p<P; p++){
        w[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (3D measure performances)
template<class T>
interp<T>::interp(T xinp[], T yinp[], T zinp[], T winp[], T xvalinp[], T yvalinp[], T zvalinp[], int Ninp, int Minp, int Linp, int Pinp, int Qinp, int Rinp, bool meas){
    D  = 3;
    measure = meas;
    N = Ninp;
    M = Minp;
    L = Linp;
    P = Pinp;
    Q = Qinp;
    R = Rinp;
    x.resize(N);
    y.resize(M);
    z.resize(L);
    xval.resize(P);
    yval.resize(Q);
    zval.resize(R);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int m=0; m<M; m++){
        y[m] = yinp[m];
    }
    for(unsigned int l=0; l<L; l++){
        z[l] = zinp[l];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[n][m][l] = winp[l+m*L+n*M*L]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int q=0; q<Q; q++){
        yval[q] = yvalinp[q];
    }
    for(unsigned int r=0; r<R; r++){
        zval[r] = zvalinp[r];
    }
    for(unsigned int p=0; p<P; p++){
        w[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
  
// Constructor (3D)
template<class T>
interp<T>::interp(vector<T> xinp, vector<T> yinp, vector<T> zinp, vector<vector<vector<T> > > winp, vector<T> xvalinp, vector<T> yvalinp, vector<T> zvalinp){
    D  = 3;
    N = xinp.size();
    M = yinp.size();
    L = zinp.size();
    P = xvalinp.size();
    Q = yvalinp.size();
    R = zvalinp.size();
    x.resize(N);
    y.resize(M);
    z.resize(L);
    xval.resize(P);
    yval.resize(Q);
    zval.resize(R);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int m=0; m<M; m++){
        y[m] = yinp[m];
    }
    for(unsigned int l=0; l<L; l++){
        z[l] = zinp[l];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[n][m][l] = winp[l+m*L+n*M*L]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int q=0; q<Q; q++){
        yval[q] = yvalinp[q];
    }
    for(unsigned int r=0; r<R; r++){
        zval[r] = zvalinp[r];
    }
    for(unsigned int p=0; p<P; p++){
        w[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}
    
// Constructor (3D measure performances)
template<class T>
interp<T>::interp(vector<T> xinp, vector<T> yinp, vector<T> zinp, vector<vector<vector<T> > > winp, vector<T> xvalinp, vector<T> yvalinp, vector<T> zvalinp, bool meas){
    D  = 3;
    measure = meas;
    N = xinp.size();
    M = yinp.size();
    L = zinp.size();
    P = xvalinp.size();
    Q = yvalinp.size();
    R = zvalinp.size();
    x.resize(N);
    y.resize(M);
    z.resize(L);
    xval.resize(P);
    yval.resize(Q);
    zval.resize(R);
    w.resize(K);
    f.resize(P);
    for(unsigned int n=0; n<N; n++){
        x[n] = xinp[n];
    }
    for(unsigned int m=0; m<M; m++){
        y[m] = yinp[m];
    }
    for(unsigned int l=0; l<L; l++){
        z[l] = zinp[l];
    }
    for(unsigned int k=0; k<K; k++){
        w[k].resize(N);
        for(unsigned int n=0; n<N; n++){
            w[n].resize(M);
            for(unsigned int m=0; m<M; m++){
                w[n][m].resize(L);
                for(unsigned int l=0; l<L; l++){
                    w[n][m][l] = winp[l+m*L+n*M*L]; 
                }
            }
        }
    }
    for(unsigned int p=0; p<P; p++){
        xval[p] = xvalinp[p];
    }
    for(unsigned int q=0; q<Q; q++){
        yval[q] = yvalinp[q];
    }
    for(unsigned int r=0; r<R; r++){
        zval[r] = zvalinp[r];
    }
    for(unsigned int p=0; p<P; p++){
        w[p].resize(Q);
        for(unsigned int q=0; q<Q; q++){
            f[p][q].resize(R);
            for(unsigned int r=0; r<R; r++){
                f[p][q][r] = 0; 
            }
        }
    }
}

//-----------------
//   Destructor 
//-----------------

template<class T>
interp<T>::~interp(){}

//-----------------
// General tools 
//-----------------

template<class T>
vector<vector<T> > interp<T>::NewtonDividedDifferences(vector<T> xinp, vector<T> yinp){
    vector<vector<T> > d;
    int m = yinp.size();
    d.resize(m);
    for(unsigned int i=0; i<m; i++){
        d[i].resize(m);
        d[i][0] = yinp[i];
        for(unsigned int j=1; j<m; j++){
            d[i][j] = 0;
        }
    }
    for(unsigned int j=1; j<m; j++){
        for(unsigned int i=j; i<m; i++){
            d[i][j] = (d[i][j-1]-d[i-1][j-1])/(x[i]-x[i-j+1]);
        }
    }
    return d;
}

template<class T>
T interp<T>::LagrangeInterpolant(vector<T> xinp, T xvalinp, unsigned int index){
    int n = xinp.size();
    int m = 1;
    T I;
    I = 1;
    for(unsigned int i=0; i<n; i++){
        if(i!=index){
            I *= (xvalinp-xinp[i])/(xinp[index]-xinp[i]);
        }
    }
    return I;
}

template<class T>
vector<T> interp<T>::LagrangeInterpolant(vector<T> xinp, vector<T> xvalinp, unsigned int index){
    int n = xinp.size();
    int m = xvalinp.size();
    vector<T> I;
    I.resize(m);
    for(unsigned int k=0; k<m; k++){
        I[k] = 1;
        for(unsigned int i=0; i<n; i++){
            if(i!=index){
                I[k] *= (xvalinp[k]-xinp[i])/(xinp[index]-xinp[i]);
            }
        }
    }
    return I;
}

template<class T>
void interp<T>::Lagrange1D(){
    vector<T> yinp, omega;
    yinp.resize(N);
    for(unsigned int n=0; n<N; n++){
        yinp[n] = w[n][0][0];
    }
    vector<vector<T> > d = NewtonDividedDifferences(x,yinp);
    for(unsigned int n=0; n<N; n++){
        omega.clear();
        omega.resize(P);
        for(unsigned int p=0; p<P; p++){
            omega[p] = 1;
        }
        for(unsigned int j=0; j<n; j++){
            for(unsigned int p=0; p<P; p++){
                omega[p] *= (xval[p]-x[j]);
            }
        }
        for(unsigned int p=0; p<P; p++){
            f[p][0][0] += omega[p]*d[n][n]; 
        }
    }
}

template<class T>
void interp<T>::Lagrange2D(){
    for(unsigned int q=0; q<Q; q++){
        for(unsigned int p=0; p<P; p++){
            for(unsigned int n=0; n<N; n++){
                for(unsigned int m=0; m<M; m++){
                    f[p][q][0] += w[n][m][0]*LagrangeInterpolant(x,xval[p],n)*LagrangeInterpolant(y,yval[q],m);
                }
            }
        }
    }
}

template<class T>
void interp<T>::Lagrange3D(){
    for(unsigned int r=0; r<R; r++){
        for(unsigned int q=0; q<Q; q++){
            for(unsigned int p=0; p<P; p++){
                for(unsigned int l=0; l<L; l++){
                    for(unsigned int m=0; m<M; m++){
                        for(unsigned int n=0; n<N; n++){
                            f[p][q][r] += w[n][m][l]*LagrangeInterpolant(x,xval[p],n)*LagrangeInterpolant(y,yval[q],m)*LagrangeInterpolant(z,zval[r],l);
                        }
                    }
                }
            }
        }
    }
}

template<class T>
void interp<T>::Hermite1D(){
    
}

template<class T>
void interp<T>::Hermite2D(){
    
}

template<class T>
void interp<T>::Hermite3D(){
    
}

template<class T>
void interp<T>::interpolation(int type){
    switch(type){
        case 0:
            switch(D){
                case 1:
                    Lagrange1D();
                    break;
                case 2:
                    Lagrange2D();
                    break;
                case 3:
                    Lagrange3D();
                    break;
                default:
                    cout << "Error: wrong dimension. Process interrupted." << endl;
                    break;
                
            }
            break;
        case 1:
            switch(D){
                case 1:
                    Hermite1D();
                    break;
                case 2:
                    Hermite2D();
                    break;
                case 3:
                    Hermite3D();
                    break;
                default:
                    cout << "Error: wrong dimension. Process interrupted." << endl;
                    break;
                
            }
            break;
        default:
            cout << "Error: type not recognised. Process interrupted." << endl;
            break;
    }
}

//-----------------
// API functions 
//-----------------

template<class T>
vector<T> interp<T>::f1D(){
    vector<T> res;
    if(D==1){
        res.resize(P);
        for(unsigned int p=0; p<P; p++){
            res[p] = f[p][0][0];
        }
    }else{
        cout << "Error: the interpolated functon is not 1D. Result left to zero." << endl;
    }
    return res;
}

template<class T>
vector<vector<T> > interp<T>::f2D(){
    vector<vector<T> > res;
    if(D==2){
        res.resize(P);
        for(unsigned int p=0; p<P; p++){
            res.resize(Q);
            for(unsigned int q=0; q<Q; q++){
                res[p][q] = f[p][q][0];
            }
        }
    }else{
        cout << "Error: the interpolated functon is not 2D. Result left to zero." << endl;
    }
    return res;
}

template<class T>
vector<vector<vector<T> > > interp<T>::f3D(){
    vector<vector<vector<T> > > res;
    if(D==3){
        res = f;
    }else{
        cout << "Error: the interpolated functon is not 3D. Result left to zero." << endl;
    }
    return res;
}
    
//--------------------------
// Overloading of operators 
//--------------------------

template<class U>
ostream& operator<<(ostream &output, const interp<U>& rhs){
    switch(rhs.D){
        case 1:
            for(unsigned int p=0; p<rhs.P; p++){
                output << rhs.f[p][0][0];
                if(p==rhs.P-1){
                    output << endl;
                }else{
                    output << ", ";
                }
            }
            break;
        case 2:
            for(unsigned int q=0; q<rhs.Q; q++){
                for(unsigned int p=0; p<rhs.P; p++){
                    output << rhs.f[p][q][0];
                    if(p==rhs.P-1){
                        output << endl;
                    }else{
                        output << ", ";
                    }
                }
            }
            break;
        case 3:
            for(unsigned int r=0; r<rhs.R; r++){
                for(unsigned int q=0; q<rhs.Q; q++){
                    for(unsigned int p=0; p<rhs.P; p++){
                        output << rhs.f[p][q][r];
                        if(p==rhs.P-1){
                            output << endl;
                        }else if(p==rhs.P-1 && q=rhs.Q-1){
                            output << endl;
                            output << "//" << endl;
                        }else{
                            output << ", ";
                        }
                    }
                }
            }
            break;
        default:
            cout << "Error: dimension not recognised." << endl;
            break;
    }
}