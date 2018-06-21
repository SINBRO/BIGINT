//
// Created by andrey on 20.06.18.
//

#ifndef BIGINT_VECTOR_H
#define BIGINT_VECTOR_H

#include <iostream>
#include <memory>

template<class T>
struct vector<T> {
    vector<T>();

    ~vector<T>();

    T operator[(uint32_t i);


private:
    union data {
        uint32_t * small_data[2];
        std::unique_ptr<vector<uint32_t>> dynamic_data;
    };

};

#endif //BIGINT_VECTOR_H
