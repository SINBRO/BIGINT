//
// Created by PlatinumGod on 03.05.2018.
//
#include <iostream>
#include "big_integer.h"
#include <algorithm>

big_integer::big_integer()
        : negative(false),
          data{0} {
}

big_integer::big_integer(big_integer const &other) = default;

big_integer::big_integer(int32_t a)
        : negative(a < 0),
          data{a < 0 ? -(uint32_t) a : (uint32_t) a} {
}

big_integer::big_integer(std::string const &str) {
    *this = big_integer();
    uint32_t pos = 0;
    if (str[0] == '-') {
        negative = true;
        ++pos;
    }
    for (; pos < str.size(); ++pos) {
        *this *= 10;
        add_u_int_to_abs((uint32_t) (str[pos] - '0'));
    }
    if (data.size() == 1 && data[0] == 0) {
        negative = false;
    }
}

big_integer::~big_integer() = default;

big_integer &big_integer::operator=(big_integer const &other) = default;

big_integer &big_integer::operator=(int32_t value) {
    *this = big_integer(value);
    return *this;
}

template<class FunctorT1, class FunctorT2, class FunctorT3>
inline void big_integer::sum_u_int_with_abs_at_pos(uint32_t rhs, uint32_t pos, FunctorT1 operation, FunctorT2 trigger,
                                                   uint32_t check_value, FunctorT3 operation_if_carry) {
    bool carry;
    int64_t sum = operation(data[pos], rhs);
    uint32_t ins_value = UINT32_MAX - check_value;
    if (trigger(sum)) {
        carry = true;
        sum = operation(sum, -((int64_t) UINT32_MAX) - 1);
        for (uint32_t i = pos + 1; i < data.size(); ++i) {
            if (data[i] == check_value) {
                data[i] = ins_value;
            } else {
                data[i] = operation(data[i], 1);
                carry = false;
                break;
            }
        }
        if (carry) {
            operation_if_carry(data);
        }
    }
    data[pos] = (uint32_t) sum;
}

void big_integer::add_u_int_to_abs(uint32_t rhs) {
    sum_u_int_with_abs_at_pos(rhs, 0, [](int64_t a, int64_t b) { return a + b; },
                              [](int64_t x) { return x > UINT32_MAX; }, UINT32_MAX,
                              [](std::vector<uint32_t> &arr) { arr.push_back(1); });
}

void big_integer::sub_u_int_from_abs(uint32_t rhs) {
    if (data.size() == 1 && data[0] < rhs) {
        data[0] = rhs - data[0];
        negative ^= true;
    } else {
        sum_u_int_with_abs_at_pos(rhs, 0, [](int64_t a, int64_t b) { return a - b; },
                                  [](int64_t x) { return x < 0; }, 0, [](std::vector<uint32_t> arr) {});
    }
}

big_integer &big_integer::plus_other_sign(big_integer const &rhs) {
    if (abs_bigger(*this, rhs)) {
        sub_less(rhs);
    } else {
        data = big_integer(rhs).sub_less(*this).data;
        if (*this == (uint32_t) 0) {
            negative = false;
        } else {
            negative ^= true;
        }
    }
    return *this;
}

big_integer &big_integer::plus_same_sign(big_integer const &rhs) {
    add_or_sub_less(rhs, [](int64_t a, int64_t b) { return a + b; }, [](int64_t x) { return x > UINT32_MAX; },
                    UINT32_MAX,
                    [](std::vector<uint32_t> arr) { arr.push_back(1); });
    return *this;
}

template<class FunctorT1, class FunctorT2, class FunctorT3>
inline void big_integer::add_or_sub_less(big_integer const &rhs, FunctorT1 operation, FunctorT2 trigger,
                                         uint32_t check_value, FunctorT3 operation_if_carry) {
    bool carry = false;
    int64_t inter_res;
    expand_to_size((uint32_t) rhs.data.size() + 1);
    for (uint32_t i = 0; i < rhs.data.size(); ++i) {
        inter_res = operation(data[i], (int64_t) rhs.data[i] + carry);
        if (trigger(inter_res)) {
            carry = true;
            inter_res = operation(inter_res, -(int64_t) UINT32_MAX - 1);
        } else {
            carry = false;
        }
        data[i] = (uint32_t) inter_res;
    }
    if (carry) {
        sum_u_int_with_abs_at_pos(1, (uint32_t) rhs.data.size(), operation, trigger, check_value, operation_if_carry);
    }
    delete_zero();
}

big_integer &big_integer::sub_less(big_integer const &rhs) {
    add_or_sub_less(rhs, [](int64_t a, int64_t b) { return a - b; }, [](int64_t x) { return x < 0; }, 0,
                    [](std::vector<uint32_t> arr) {});
    return *this;
}

big_integer &big_integer::operator+=(big_integer const &rhs) {
    if (negative == rhs.negative) {
        plus_same_sign(rhs);
    } else {
        plus_other_sign(rhs);
    }
    return *this;
}

big_integer &big_integer::operator+=(uint32_t rhs) {
    if (negative) {
        sub_u_int_from_abs(rhs);
    } else {
        add_u_int_to_abs(rhs);
    }
    return *this;
}

big_integer &big_integer::operator+=(int32_t rhs) {
    if (rhs < 0) {
        *this -= (uint32_t) -rhs;
    } else {
        *this += (uint32_t) rhs;
    }
    return *this;
}

big_integer &big_integer::operator-=(uint32_t rhs) {
    if (negative) {
        add_u_int_to_abs(rhs);
    } else {
        sub_u_int_from_abs(rhs);
    }
    return *this;
}

big_integer &big_integer::operator-=(int32_t rhs) {
    if (rhs < 0) {
        *this += (uint32_t) -rhs;
    } else {
        *this -= (uint32_t) rhs;
    }
    return *this;
}

big_integer &big_integer::operator-=(big_integer const &rhs) {
    if (negative != rhs.negative) {
        plus_same_sign(rhs);
    } else {
        plus_other_sign(rhs);
    }
    return *this;
}

big_integer &big_integer::operator*=(int32_t rhs) {
    negative ^= rhs < 0;
    *this *= rhs < 0 ? - (uint32_t) rhs : (uint32_t) rhs;
    return *this;
}

big_integer &big_integer::operator*=(uint32_t rhs) {
    if (rhs == 0) {
        return *this = 0;
    }
    data.push_back(0);
    std::vector<uint32_t> data_copy = data;
    data = {0};

    for (uint32_t j = 0; j < data_copy.size(); ++j) {
        data.push_back(0);
    }
    uint64_t mul = 0, sum;
    for (uint32_t i = 0; i < data_copy.size(); ++i) {
        sum = mul;
        mul = (uint64_t) data_copy[i] * rhs;
        sum += ((uint64_t) data[i]) + (mul & (uint64_t) UINT32_MAX);
        mul >>= 32;
        if (sum > UINT32_MAX) {
            mul += sum / ((uint64_t) UINT32_MAX + 1);
            sum %= ((uint64_t) UINT32_MAX + 1);
        }
        data[i] = (uint32_t) sum;
    }
    data[data.size() - 1] = (uint32_t) mul;
    delete_zero();
    return *this;
}

big_integer &big_integer::operator*=(big_integer const &rhs) {
    big_integer res = big_integer();
    res.negative = negative ^ rhs.negative;
    uint64_t mul = 0, sum;
    data.push_back(0);
    for (uint32_t k = 0; k < rhs.data.size() + data.size(); ++k) {
        res.data.push_back(0);
    }
    for (uint32_t i = 0; i < rhs.data.size(); ++i) {
        for (uint32_t j = 0; j < data.size(); ++j) {
            sum = mul;
            mul = (uint64_t) data[j] * rhs.data[i];
            sum += (uint64_t) res.data[i + j] + (mul & (uint64_t) UINT32_MAX);
            mul >>= 32;
            if (sum > UINT32_MAX) {
                mul += sum / ((uint64_t) UINT32_MAX + 1);
                sum %= ((uint64_t) UINT32_MAX + 1);
            }
            res.data[i + j] = (uint32_t) sum;
        }
    }
    delete_zero();
    res.delete_zero();
    *this = res;
    return *this;
}

big_integer &big_integer::operator/=(uint32_t rhs) {
    return *this = divide(rhs);
}

big_integer &big_integer::operator/=(int32_t rhs) {
    if (rhs < 0) {
        negative ^= true;
        rhs = -rhs;
    }
    return *this /= (uint32_t) rhs;
}

big_integer &big_integer::operator%=(uint32_t rhs) {
    divide(rhs);
    return *this;
}

inline big_integer big_integer::divide(uint32_t rhs) {
    if (data.size() == 1 && data[0] < rhs) {
        return big_integer();
    }

    big_integer res = big_integer();
    res.expand_to_size((uint32_t) data.size());
    res.negative = negative;

    uint64_t inter_quotient = 0;
    uint32_t inter_dividend;

    for (int32_t i = (uint32_t) data.size() - 1; i >= 0; --i) {
        inter_quotient = (inter_quotient << 32) + data[i];
        inter_dividend = (uint32_t) (inter_quotient / rhs);
        inter_quotient %= rhs;
        res.data[i] = inter_dividend;
    }

    res.delete_zero();
    delete_zero();
    data = {(uint32_t) inter_quotient};
    return res;
}

big_integer &big_integer::operator/=(big_integer const &rhs) {
    if (rhs.data.size() == 1) {
        negative ^= rhs.negative;
        return *this /= rhs.data[0];
    }
    if (abs_bigger(rhs, *this)) {
        return *this = big_integer();
    }
    *this = divide(rhs);
    return *this;
}

inline big_integer big_integer::divide(big_integer const &rhs) { //*this -> *this % rhs;  returns *this / rhs;
    big_integer res = big_integer();
    big_integer rhs_copy = rhs;
    big_integer inter_part;

    auto normalization_factor = (uint32_t) (((int64_t) UINT32_MAX + 1) / ((int64_t) rhs.data.back() + 1));
    rhs_copy *= normalization_factor;
    *this *= normalization_factor;

    data.push_back(0);

    res.expand_to_size((uint32_t) (data.size() - rhs_copy.data.size() + 1));
    res.negative = negative ^ rhs.negative;

    uint64_t inter_dividend;
    uint32_t inter_quotient;
    uint32_t inter_divider = rhs_copy.data.back();

    for (uint32_t pos = (uint32_t) data.size() - 1; pos >= rhs_copy.data.size(); --pos) {
        inter_dividend = ((uint64_t) data[pos] << 32) + data[pos - 1];
        inter_quotient = (uint32_t) (inter_dividend /
                                     inter_divider); // if it can be < (feel like no) than it should be....

        inter_part = rhs_copy * inter_quotient;
        inter_part.expand_to_size((uint32_t) rhs_copy.data.size() + 1);
        while (compare_at_pos(inter_part, pos) == -1) {
            --inter_quotient;
            inter_part -= rhs_copy;
        }
        sub_at_pos(inter_part, pos);
        res.data[pos - rhs_copy.data.size()] = inter_quotient;
    }
    res.delete_zero();
    *this /= normalization_factor;
    if (data == std::vector<uint32_t>({0})) {
        negative = false;
    }
    return res;
}

inline void big_integer::sub_at_pos(big_integer const &rhs, int32_t pos) { // no carry after last digit
    bool carry = false;
    int64_t sub;
    pos -= rhs.data.size() - 1;
    for (uint32_t i : rhs.data) {
        sub = (int64_t) data[pos] - i - carry;
        if (sub < 0) {
            carry = true;
            sub = sub + UINT32_MAX + 1;
        } else {
            carry = false;
        }
        data[pos] = (uint32_t) sub;
        ++pos;
    }
}

inline int32_t big_integer::compare_at_pos(big_integer const &rhs, int32_t pos) const {
    pos++;
    for (uint32_t i = 1; i <= rhs.data.size(); ++i) {
        if (data[pos - i] < rhs.data[rhs.data.size() - i]) {
            return -1;
        }
        if (data[pos - i] > rhs.data[rhs.data.size() - i]) {
            return 1;
        }
    }
    return 0;
}

big_integer &big_integer::operator%=(big_integer const &rhs) {
    /*if (rhs == 1 || rhs == -1) {
        return *this = 0;
    }*/
    if (abs_bigger(rhs, *this)) {
        return *this;
    }
    divide(rhs);
    return *this;
}

inline void big_integer::delete_zero() {
    while (data.back() == 0 && data.size() > 1) {
        data.pop_back();
    }
}

inline void big_integer::expand_to_size(uint32_t new_size) {
    while (data.size() < new_size) {
        data.push_back(0);
    }
}

void big_integer::change_form() {
    if (negative) {
        bool trigger = false;
        for (uint32_t &d: data) {
            if (!trigger) {
                d = (uint32_t) -d;
                if (d != 0) {
                    trigger = true;
                }
            } else {
                d = ~d;
            }
        }
    }
}

template<class FunctorT>
big_integer &big_integer::apply_bitwise_operation(big_integer const &rhs, FunctorT functor) {
    big_integer rhs_copy = rhs;

    data.push_back(0);
    expand_to_size((uint32_t) rhs_copy.data.size());
    rhs_copy.expand_to_size((uint32_t) data.size());

    change_form();
    rhs_copy.change_form();

    for (uint32_t i = 0; i < data.size(); ++i) {
        data[i] = functor(data[i], rhs_copy.data[i]);
    }

    if ((data.back() & (1 << 31)) != 0) {
        negative = true;
    }
    change_form();
    delete_zero();

    return *this;
}

big_integer &big_integer::operator&=(big_integer const &rhs) {
    return apply_bitwise_operation(rhs, std::bit_and<uint32_t>());
}

big_integer &big_integer::operator|=(big_integer const &rhs) {
    return apply_bitwise_operation(rhs, std::bit_or<uint32_t>());
}

big_integer &big_integer::operator^=(big_integer const &rhs) {
    return apply_bitwise_operation(rhs, std::bit_xor<uint32_t>());
}

template<class FunctorT1, class FunctorT2, class FunctorT3>
inline void big_integer::shift(uint32_t rhs, FunctorT1 shift_big_digits, FunctorT2 shift_digits, FunctorT3 final_action) {
    if (*this == 0) {
        return;
    }
    uint32_t fulls = rhs / 32;
    rhs %= 32;

    auto least_data = shift_big_digits(data, fulls);

    shift_digits(least_data, rhs);

    data.insert(data.end(), least_data.begin(), least_data.end());

    if (data.empty()) {
        data.push_back(0);
    }
    delete_zero();
    final_action(*this);
}

big_integer &big_integer::operator<<=(uint32_t rhs) {
    shift(rhs,
          [](std::vector<uint32_t> &data, uint32_t fulls) {
              auto res = data;
              data.clear();
              for (uint32_t i = 0; i < fulls; ++i) {
                  data.push_back(0);
              }
              return res;
          },
          [](std::vector<uint32_t> &data, uint32_t rhs) {
              uint64_t shifted = 0;
              data.push_back(0);
              for (uint32_t &i : data) {
                  shifted = (((uint64_t) i) << rhs) + (shifted >> 32);
                  i = (uint32_t) shifted;
              }
              if (data.back() == 0) {
                  data.pop_back();
              }
          },
          [](big_integer x) {});
    return *this;
}

big_integer &big_integer::operator>>=(uint32_t rhs) {
    shift(rhs,
          [](std::vector<uint32_t> &data, uint32_t fulls) {
              auto res = std::vector<uint32_t>();
              for (uint32_t j = fulls; j < data.size(); ++j) {
                  res.push_back(data[j]);
              }
              data.clear();
              return res;
          },
          [](std::vector<uint32_t> &data, uint32_t rhs) {
              uint32_t rhs_mask = ~(UINT32_MAX << rhs);
              uint32_t part = 0, new_digit;
              for (int32_t i = (uint32_t) data.size() - 1; i >= 0; --i) {
                  new_digit = (part << (32 - rhs)) + (data[i] >> rhs);
                  part = data[i] & rhs_mask;
                  data[i] = new_digit;
              }
          },
          [](big_integer &x) {
              if (x.negative) {
                  --x;
              }
          });
    return *this;
}

big_integer &big_integer::operator>>=(int32_t rhs) {
    if (rhs < 0) {
        *this <<= (uint32_t) -rhs;
    } else {
        *this >>= (uint32_t) rhs;
    }
    return *this;
}

big_integer &big_integer::operator<<=(int32_t rhs) {
    if (rhs < 0) {
        *this >>= (uint32_t) -rhs;
    } else {
        *this <<= (uint32_t) rhs;
    }
    return *this;
}

big_integer big_integer::operator+() const {
    return *this;
}

big_integer big_integer::operator-() const {
    big_integer res = *this;
    if (res != 0) {
        res.negative ^= true;
    }
    return res;
}

big_integer big_integer::operator~() const {
    return -(*this + (uint32_t) 1);
}

big_integer &big_integer::operator++() {
    return *this += 1;
}

const big_integer big_integer::operator++(int) {
    big_integer r = *this;
    *this += 1;
    return r;
}

big_integer &big_integer::operator--() {
    return *this -= 1;
}

const big_integer big_integer::operator--(int) {
    big_integer r = *this;
    *this -= 1;
    return r;
}

big_integer operator+(big_integer a, big_integer const &b) {
    return a += b;
}

big_integer operator+(big_integer a, int32_t b) {
    return a += b;
}

big_integer operator+(big_integer a, uint32_t b) {
    return a += b;
}

big_integer operator+(int32_t a, big_integer b) {
    return b += a;
}

big_integer operator+(uint32_t a, big_integer b) {
    return b += a;
}

big_integer operator-(big_integer a, big_integer const &b) {
    return a -= b;
}

big_integer operator-(big_integer a, int32_t b) {
    return a -= b;
}

big_integer operator*(big_integer a, int32_t b) {
    return a *= b;
}

big_integer operator*(big_integer a, uint32_t b) {
    return a *= b;
}

big_integer operator*(big_integer a, big_integer const &b) {
    return a *= b;
}

big_integer operator/(big_integer a, big_integer const &b) {
    return a /= b;
}

big_integer operator%(big_integer a, big_integer const &b) {
    return a %= b;
}

big_integer operator&(big_integer a, big_integer const &b) {
    return a &= b;
}

big_integer operator|(big_integer a, big_integer const &b) {
    return a |= b;
}

big_integer operator^(big_integer a, big_integer const &b) {
    return a ^= b;
}

big_integer operator<<(big_integer a, int32_t b) {
    return a <<= b;
}

big_integer operator>>(big_integer a, int32_t b) {
    return a >>= b;
}

bool operator==(big_integer const &a, uint32_t b) {
    return a.data.size() == 1 && !a.negative && a.data[0] == b;
}

bool operator==(big_integer const &a, int32_t b) {
    return a.data.size() == 1 && a.negative == (b < 0) && a.data[0] == (b < 0 ? -(uint32_t) b : (uint32_t) b);
}

bool operator!=(big_integer const &a, int32_t b) {
    return !(a == b);
}

bool operator!=(big_integer const &a, uint32_t b) {
    return !(a == b);
}

bool operator==(big_integer const &a, big_integer const &b) {
    if (a.negative != b.negative || a.data.size() != b.data.size()) {
        return false;
    }
    for (uint32_t i = 0; i < a.data.size(); ++i) {
        if (a.data[i] != b.data[i]) {
            return false;
        }
    }
    return true;
}

bool operator!=(big_integer const &a, big_integer const &b) {
    return !(a == b);
}

bool abs_bigger(big_integer const &a, big_integer const &b) {
    if (a.data.size() != b.data.size()) {
        return a.data.size() > b.data.size();
    }
    return a.compare_at_pos(b, (uint32_t) a.data.size() - 1) == 1;
}

bool operator<(big_integer const &a, big_integer const &b) {
    if (a.negative) {
        if (!b.negative) {
            return true;
        }
        return abs_bigger(a, b);
    }
    if (b.negative) {
        return false;
    }
    return abs_bigger(b, a);
}

bool operator>(big_integer const &a, big_integer const &b) {
    return b < a;
}

bool operator<=(big_integer const &a, big_integer const &b) {
    return !(a > b);
}

bool operator>=(big_integer const &a, big_integer const &b) {
    return !(a < b);
}

std::string to_string(big_integer const &a) {
    std::vector<char> ans;
    big_integer a_copy = a;
    uint32_t d = 10;
    a_copy.negative = false;
    do {
        auto inter_res = a_copy.divide(d);
        ans.push_back((char) a_copy.data[0] + '0');
        a_copy = inter_res;
    } while (a_copy != 0);
    if (a.negative) {
        ans.push_back('-');
    }
    reverse(ans.begin(), ans.end());
    return std::string(&ans[0], ans.size());
}

std::ostream &operator<<(std::ostream &s, big_integer const &a) {
    return s << to_string(a);
}
