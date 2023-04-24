#ifndef COMPLEX_FMT_HH
#define COMPLEX_FMT_HH

#include "fmt/format.h"
#include <complex>

/// From https://gitlab.com/tesch1/cppduals/blob/master/duals/dual#L1304
/// std::complex<> Formatter for libfmt https://github.com/fmtlib/fmt
///
/// libfmt does not provide a formatter for std::complex<>, although
/// one is proposed for c++20.  Anyway, at the expense of a k or two,
/// you can define CPPDUALS_LIBFMT_COMPLEX and get this one.
///
/// The standard iostreams formatting of complex numbers is (a,b),
/// where a and b are the real and imaginary parts.  This formats a
/// complex number (a+bi) as (a+bi), offering the same formatting
/// options as the underlying type - with the addition of three
/// optional format options, only one of which may appear directly
/// after the ':' in the format spec (before any fill or align): '$'
/// (the default if no flag is specified), '*', and ','.  The '*' flag
/// adds a * before the 'i', producing (a+b*i), where a and b are the
/// formatted value_type values.  The ',' flag simply prints the real
/// and complex parts separated by a comma (same as iostreams' format).
/// As a concrete exmple, this formatter can produce either (3+5.4i)
/// or (3+5.4*i) or (3,5.4) for a complex<double> using the specs {:g}
/// | {:$g}, {:*g}, or {:,g}, respectively.  (this implementation is a
/// bit hacky - glad for cleanups).
///
template <typename T, typename Char>
struct fmt::formatter<std::complex<T>, Char> : public fmt::formatter<T, Char> {
    using base = fmt::formatter<T, Char>;
    enum style { expr, star, pair } style_ = expr;
    fmt::detail::dynamic_format_specs<Char> specs_;
    FMT_CONSTEXPR auto parse(format_parse_context &ctx) -> decltype(ctx.begin()) {
        using handler_type = fmt::detail::dynamic_specs_handler<format_parse_context>;
        auto type = fmt::detail::type_constant<T, Char>::value;
        fmt::detail::specs_checker<handler_type> handler(handler_type(specs_, ctx), type);
        auto it = ctx.begin();
        if(it != ctx.end()) {
            switch(*it) {
            case '$':
                style_ = style::expr;
                ctx.advance_to(++it);
                break;
            case '*':
                style_ = style::star;
                ctx.advance_to(++it);
                break;
            case ',':
                style_ = style::pair;
                ctx.advance_to(++it);
                break;
            default:
                break;
            }
        }
        parse_format_specs(ctx.begin(), ctx.end(), handler);
        // todo: fixup alignment
        return base::parse(ctx);
    }
    template <typename FormatCtx>
    auto format(const std::complex<T> &x, FormatCtx &ctx) -> decltype(ctx.out()) {
        format_to(ctx.out(), "(");
        if(style_ == style::pair) {
            base::format(x.real(), ctx);
            format_to(ctx.out(), ",");
            base::format(x.imag(), ctx);
            return format_to(ctx.out(), ")");
        }
        if(x.real() || !x.imag()) base::format(x.real(), ctx);
        if(x.imag()) {
            if(x.real() && x.imag() >= 0 && specs_.sign != sign::plus) format_to(ctx.out(), "+");
            base::format(x.imag(), ctx);
            if(style_ == style::star)
                format_to(ctx.out(), "*i");
            else
                format_to(ctx.out(), "i");
            if(std::is_same<typename std::decay<T>::type, float>::value) format_to(ctx.out(), "f");
            if(std::is_same<typename std::decay<T>::type, long double>::value)
                format_to(ctx.out(), "l");
        }
        return format_to(ctx.out(), ")");
    }
};

#endif
