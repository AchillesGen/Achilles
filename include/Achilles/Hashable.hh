#pragma once

#include <functional>
#include <type_traits>

namespace achilles::utils {

// From: https://stackoverflow.com/a/57595105
template <typename Head, typename... Tail>
void hash_combine(std::size_t &seed, const Head &head, const Tail &...tail) {
    std::hash<Head> hasher;
    seed ^= hasher(head) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    (hash_combine(seed, tail), ...);
}

template <typename Derived> class Hashable {
  public:
    std::size_t hash() const { return static_cast<const Derived *>(this)->hash_impl(); }
};

} // namespace achilles::utils
