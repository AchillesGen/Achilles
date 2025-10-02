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

void cmap_init(std::map<std::string, std::complex<double>> *&map) {
    map = new std::map<std::string, std::complex<double>>();
}

void cmap_insert(std::map<std::string, std::complex<double>> *map, const char *key,
                 std::complex<double> value) {
    map->insert({std::string(key), value});
}

void cmap_lookup(std::map<std::string, std::complex<double>> *map, const char *key,
                 std::complex<double> *outval) {
    std::string skey(key);
    auto it = map->find(skey);
    if(it == map->end()) {
        *outval = std::complex<double>();
    } else {
        *outval = it->second;
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
