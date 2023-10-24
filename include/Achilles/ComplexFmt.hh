#ifndef COMPLEX_FMT_HH
#define COMPLEX_FMT_HH

#include "fmt/format.h"
#include <complex>

// from https://github.com/fmtlib/fmt/issues/1467

template <typename T, typename Char>
struct fmt::formatter<std::complex<T>, Char> : public fmt::formatter<T, Char> {
  private:
    typedef fmt::formatter<T, Char> base;
    enum style { expr, star, pair } style_ = expr;

    detail::dynamic_format_specs<Char> specs_;

  public:
    template <typename ParseContext> FMT_CONSTEXPR auto parse(ParseContext &ctx) -> const Char * {
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

        auto type = detail::type_constant<T, Char>::value;
        auto end = detail::parse_format_specs(ctx.begin(), ctx.end(), specs_, ctx, type);
        if(type == detail::type::char_type) detail::check_char_specs(specs_);
        return end;
    }

    template <typename FormatContext>
    FMT_CONSTEXPR auto format(const std::complex<T> &x, FormatContext &ctx) const
        -> decltype(ctx.out()) {
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
