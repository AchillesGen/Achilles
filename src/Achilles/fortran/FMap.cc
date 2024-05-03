#include <complex>
#include <map>
#include <string>

extern "C" {

void map_init(std::map<std::string, double> *&map) {
    map = new std::map<std::string, double>();
}

void map_insert(std::map<std::string, double> *map, const char *key, double value) {
    map->insert({std::string(key), value});
}

double map_lookup(std::map<std::string, double> *map, const char *key) {
    std::string skey(key);
    //for(auto &pair : *map) { printf("key: %s, value: %f\n", pair.first.c_str(), pair.second); }
    if(map->find(skey) == map->end()) {
        return 0.0;
    } else {
        return (*map)[skey];
    }
}

bool map_contains(std::map<std::string, double> *map, const char *key) {
    std::string skey(key);
    return map->find(skey) != map->end();
}

bool map_empty(std::map<std::string, double> *map) {
    return map->empty();
}

void map_erase(std::map<std::string, double> *map, const char *key) {
    map->erase(std::string(key));
}

void map_clear(std::map<std::string, double> *map) {
    map->clear();
}

void map_delete(std::map<std::string, double> *map) {
    delete map;
}

struct complex {
    double real;
    double imag;

#ifdef __cplusplus
    complex() : real(0.0), imag(0.0) {}
    complex(double r) : real(r), imag(0.0) {}
    complex(double r, double i) : real(r), imag(i) {}
    complex(std::complex<double> c) : real(c.real()), imag(c.imag()) {}
    operator std::complex<double>() const { return std::complex<double>(real, imag); }
#endif
};

void cmap_init(std::map<std::string, std::complex<double>> *&map) {
    map = new std::map<std::string, std::complex<double>>();
}

void cmap_insert(std::map<std::string, std::complex<double>> *map, const char *key,
                 std::complex<double> value) {
    map->insert({std::string(key), value});
}

complex *cmap_lookup(std::map<std::string, std::complex<double>> *map, const char *key) {
    std::string skey(key);
    if(map->find(skey) == map->end()) {
        complex *c = new complex();
        return c;
    } else {
        return new complex((*map)[skey]);
    }
}

bool cmap_contains(std::map<std::string, std::complex<double>> *map, const char *key) {
    std::string skey(key);
    return map->find(skey) != map->end();
}

bool cmap_empty(std::map<std::string, std::complex<double>> *map) {
    return map->empty();
}

void cmap_erase(std::map<std::string, std::complex<double>> *map, const char *key) {
    map->erase(std::string(key));
}

void cmap_clear(std::map<std::string, std::complex<double>> *map) {
    map->clear();
}

void cmap_delete(std::map<std::string, std::complex<double>> *map) {
    delete map;
}
}
