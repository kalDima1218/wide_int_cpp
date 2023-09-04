#include "bits/stdc++.h"

using namespace std;

typedef complex<long double> comp;

long long nextPowerOfTwo(long long n) {
    n--;
    n |= n >> 1ll;
    n |= n >> 2ll;
    n |= n >> 4ll;
    n |= n >> 8ll;
    n |= n >> 16ll;
    n |= n >> 32ll;
    n++;
    return n;
}

void fft(vector<comp> &p, comp w, vector<vector<vector<comp>>> &buff){
    if(p.size() == 1){
        return;
    }
    long long lg = log2(p.size())-1;
    for(long long i = 0; i < p.size(); ++i){
        if(i&1){
            buff[lg][1][i/2] = p[i];
        }
        else{
            buff[lg][0][i/2] = p[i];
        }
    }
    fft(buff[lg][0], w*w, buff);
    fft(buff[lg][1], w*w, buff);
    comp wk = 1;
    long long k = p.size()/2;
    for(long long i = 0; i < p.size(); ++i){
        p[i] = buff[lg][0][i%k] + wk * buff[lg][1][i%k];
        wk*=w;
    }
}

comp get_w(long long n){
    return comp(cos(2*M_PI / (long double)n), sin(2*M_PI / (long double)n));
}

vector<comp> evalute(vector<long long> &p, vector<vector<vector<comp>>> &buff){
    vector<comp> _p(p.size());
    for(long long i = 0; i < p.size(); ++i){
        _p[i] = comp(p[i], 0);
    }
    fft(_p, get_w(p.size()), buff);
    return _p;
}

vector<long long> interpolate(vector<comp> &p, vector<vector<vector<comp>>> &buff){
    fft(p, get_w(-p.size()), buff);
    vector<long long> _p(p.size());
    for(long long i = 0; i < p.size(); ++i){
        _p[i] = round(p[i].real()/(long double)p.size());
    }
    return _p;
}

struct wide_int{
    vector<long long> val;
    bool f;
    wide_int(){}
    wide_int(long long x) {
        val = {x};
        if(x < 0) {
            f = 1;
        }
        else{
            f = 0;
        }
        carry();
    }
    long long size() {
        return val.size();
    }
    void resize(long long n) {
        while(size() < n){
            val.push_back(0);
        }
    }
    void print(){
        if(f){
            cout << "-";
        }
        long long i = val.size() - 1;
        for(;i > 0 and val[i] == 0; i--){}
        for(; i >= 0; --i){
            cout << val[i];
        }
    }
    void carry() {
        long long size = val.size();
        for (long long i = 0; i < size-1; i++){
            val[i+1] += val[i] / 10;
            val[i] %= 10;
        }
        while(abs(val[size-1]) > 9) {
            val.push_back(val[size-1]/10);
            size++;
            val[size-2] %= 10;
        }
        while(size > 1 && val[size-1] == 0) {
            val.pop_back();
            size--;
        }
        resize(nextPowerOfTwo(size));
    }
};

bool operator<(wide_int a, wide_int b) {
    if(a.f != b.f) {
        return a.f;
    }
    long long sz = nextPowerOfTwo(max(a.size(), b.size()));
    a.resize(sz);
    b.resize(sz);
    for(long long i = sz - 1; i >= 0; i--) {
        if(a.val[i] < b.val[i]) {
            return !a.f;
        }
        if(a.val[i] > b.val[i]) {
            return a.f;
        }
    }
    return a.f;
}

bool operator>(wide_int a, wide_int b) {
    if(a.f != b.f) {
        return !a.f;
    }
    long long sz = nextPowerOfTwo(max(a.size(), b.size()));
    a.resize(sz);
    b.resize(sz);
    for(long long i = sz - 1; i >= 0; i--) {
        if(a.val[i] > b.val[i]) {
            return !a.f;
        }
        if(a.val[i] < b.val[i]) {
            return a.f;
        }
    }
    return a.f;
}

bool operator==(wide_int a, wide_int b) {
    if(a.f != b.f) {
        return false;
    }
    long long sz = nextPowerOfTwo(max(a.size(), b.size()));
    a.resize(sz);
    b.resize(sz);
    for(long long i = sz - 1; i >= 0; i--) {
        if(a.val[i] != b.val[i]) {
            return false;
        }
    }
    return true;
}

bool operator<=(wide_int a, wide_int b) {
    return a < b or a == b;
}

bool operator>=(wide_int a, wide_int b) {
    return a > b or a == b;
}

wide_int operator+(wide_int a, wide_int b) {
    long long sz = nextPowerOfTwo(max(a.size(), b.size()) + 1);
    a.resize(sz);
    b.resize(sz);
    wide_int res = 0;
    res.resize(sz);
    for(long long i = 0; i < sz; i++) {
        res.val[i] = (1-a.f*2)*a.val[i] + (1-b.f*2)*b.val[i];
    }
    res.carry();
    sz = res.size();
    for(long long i = 0; i < sz-1; i++) {
        if(res.val[i] < 0) {
            res.val[i] += 10;
            res.val[i+1]--;
        }
    }
    for(long long i = 0; i < sz; i++) {
        if(res.val[i] < 0) {
            res.val[i] = -res.val[i];
            res.f = (res.f + 1) % 2;
            break;
        }
    }
    return res;
}

wide_int operator-(wide_int a, wide_int b) {
    long long sz = nextPowerOfTwo(max(a.size(), b.size()) + 1);
    a.resize(sz);
    b.resize(sz);
    wide_int res = 0;
    res.resize(sz);
    for(long long i = 0; i < sz; i++) {
        res.val[i] = (1-a.f*2)*a.val[i] - (1-b.f*2)*b.val[i];
    }
    res.carry();
    sz = res.size();
    for(long long i = 0; i < sz-1; i++) {
        if(res.val[i] < 0) {
            res.val[i] += 10;
            res.val[i+1]--;
        }
    }
    for(long long i = 0; i < sz; i++) {
        if(res.val[i] < 0) {
            res.val[i] = -res.val[i];
            res.f = (res.f + 1) % 2;
            break;
        }
    }
    return res;
}

wide_int operator*(wide_int a, wide_int b){
    long long sz = nextPowerOfTwo(a.size() + b.size() + 1);
    a.resize(sz);
    b.resize(sz);
    vector<vector<vector<comp>>> buff((long long)(log2(sz))+1, vector<vector<comp>>(2));
    for(long long i = 0; i < buff.size(); ++i){
        buff[i][0].resize(1ll<<i);
        buff[i][1].resize(1ll<<i);
    }
    vector<comp> _a = evalute(a.val, buff);
    vector<comp> _b = evalute(b.val, buff);
    vector<comp> _res(a.val.size());
    for(long long i = 0; i < _res.size(); ++i){
        _res[i] = _a[i] * _b[i];
    }
    vector<long long> res = interpolate(_res, buff);
    wide_int wi_res;
    wi_res.val = res;
    wi_res.f = a.f ^ b.f;
    wi_res.carry();
    return wi_res;
}

wide_int operator/(wide_int a, wide_int b){
    long long base = 2;
    wide_int res = 0;
    bool f = a.f ^ b.f;
    a.f = 0;
    b.f = 0;
    while(b <= a) {
        wide_int _b = b, tmp = 1;
        while(_b * base <= a){
            _b = _b * base;
            tmp = tmp * base;
        }
        a = a - _b;
        res = res + tmp;
    }
    res.f = f;
    return res;
}

wide_int operator%(wide_int a, wide_int b){
    return a - (b * (a / b));
}

wide_int operator^(wide_int a, wide_int n){
    wide_int res = 1;
    while(n > 0){
        if(n % 2 == 1){
            res = res * a;
        }
        a = a * a;
        n = n / 2;
    }
    return res;
}

signed main(){
    wide_int a, b, c;
    a = wide_int(2);
    b = wide_int(128);
    c = a^b;
    c.print();
}
