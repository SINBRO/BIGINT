//
// Created by PlatinumGod on 03.05.2018.
//

#ifndef BIGINTEGER_BIGINTEGER_H
#define BIGINTEGER_BIGINTEGER_H

#include <iostream>
#include <vector>

struct big_integer {
    big_integer();

    big_integer(big_integer const &other);

    big_integer(int32_t a);

    explicit big_integer(std::string const &str);

    ~big_integer();

    big_integer &operator=(big_integer const &other);

    big_integer &operator=(int32_t value);

    big_integer &operator+=(big_integer const &rhs);

    big_integer &operator+=(uint32_t rhs);

    big_integer &operator+=(int32_t rhs);

    big_integer &operator-=(big_integer const &rhs);

    big_integer &operator-=(uint32_t rhs);

    big_integer &operator-=(int32_t rhs);

    big_integer &operator*=(big_integer const &rhs);

    big_integer &operator*=(int32_t rhs);

    big_integer &operator*=(uint32_t rhs);

    big_integer &operator/=(big_integer const &rhs);

    big_integer &operator/=(uint32_t rhs);

    big_integer &operator%=(uint32_t rhs);

    big_integer &operator/=(int32_t rhs);

    big_integer &operator%=(big_integer const &rhs);

    big_integer &operator&=(big_integer const &rhs);

    big_integer &operator|=(big_integer const &rhs);

    big_integer &operator^=(big_integer const &rhs);

    big_integer &operator<<=(int32_t rhs);

    big_integer &operator>>=(int32_t rhs);

    big_integer &operator<<=(uint32_t rhs);

    big_integer &operator>>=(uint32_t rhs);

    big_integer operator+() const;

    big_integer operator-() const;

    big_integer operator~() const;

    big_integer &operator++();

    const big_integer operator++(int);

    big_integer &operator--();

    const big_integer operator--(int);

    friend bool operator==(big_integer const &a, big_integer const &b);

    friend bool operator!=(big_integer const &a, big_integer const &b);

    friend bool operator<(big_integer const &a, big_integer const &b);

    friend bool operator>(big_integer const &a, big_integer const &b);

    friend bool operator<=(big_integer const &a, big_integer const &b);

    friend bool operator>=(big_integer const &a, big_integer const &b);

    friend std::string to_string(big_integer const &a);

    friend bool operator==(big_integer const &a, uint32_t b);

    friend bool operator==(big_integer const &a, int32_t b);

    friend bool operator!=(big_integer const &a, uint32_t b);

    friend bool operator!=(big_integer const &a, int32_t b);

    friend big_integer operator*(big_integer a, int32_t b);

    friend big_integer operator*(big_integer a, uint32_t b);

private:
    bool negative;

    std::vector<uint32_t> data;

    big_integer &plus_same_sign(big_integer const &rhs);

    big_integer &plus_other_sign(big_integer const &rhs);

    big_integer &sub_less(big_integer const &rhs);

    void delete_zero();

    friend bool abs_bigger(big_integer const &a, big_integer const &b);

    big_integer divide(big_integer const &rhs);

    big_integer divide(uint32_t);

    void expand_to_size(uint32_t new_size);

    void sub_at_pos(big_integer const &rhs, int32_t pos);

    int32_t compare_at_pos(big_integer const &rhs, int32_t pos) const; // 1 if *this is greater

    template<class FunctorT1, class FunctorT2, class FunctorT3>
    void sum_u_int_with_abs_at_pos(uint32_t rhs, uint32_t pos, FunctorT1 operation, FunctorT2 trigger,
                                   uint32_t check_value,
                                   FunctorT3 operation_if_carry);

    template<class FunctorT1, class FunctorT2, class FunctorT3>
    void add_or_sub_less(big_integer const &rhs, FunctorT1 operation, FunctorT2 trigger, uint32_t check_value,
                         FunctorT3 operation_if_carry);

    void add_u_int_to_abs(uint32_t rhs);

    void sub_u_int_from_abs(uint32_t rhs);

    void change_form();

    template<class FunctorT>
    big_integer &apply_bitwise_operation(big_integer const &rhs, FunctorT functor);

    template<class FunctorT1, class FunctorT2, class FunctorT3>
    void shift(uint32_t rhs, FunctorT1 shift_big_digits, FunctorT2 shift_digits, FunctorT3 final_action);
};

big_integer operator+(big_integer a, big_integer const &b);

big_integer operator+(big_integer a, int32_t b);

big_integer operator+(int32_t a, big_integer b);

big_integer operator+(big_integer a, uint32_t b);

big_integer operator+(uint32_t a, big_integer b);

big_integer operator-(big_integer a, big_integer const &b);

big_integer operator-(big_integer a, int32_t b);

big_integer operator*(big_integer a, big_integer const &b);

big_integer operator/(big_integer a, big_integer const &b);

big_integer operator%(big_integer a, big_integer const &b);

big_integer operator&(big_integer a, big_integer const &b);

big_integer operator|(big_integer a, big_integer const &b);

big_integer operator^(big_integer a, big_integer const &b);

big_integer operator<<(big_integer a, int32_t b);

big_integer operator>>(big_integer a, int32_t b);

big_integer operator<<(big_integer a, uint32_t b);

big_integer operator>>(big_integer a, uint32_t b);

std::ostream &operator<<(std::ostream &s, big_integer const &a);

#endif //BIGINTEGER_BIGINTEGER_H
