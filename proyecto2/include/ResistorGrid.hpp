//
// Created by ubuntu on 01/10/18.
//

#ifndef PROYECTO2_RESISTORGRID_H
#define PROYECTO2_RESISTORGRID_H


namespace anpi{
  /// Pack a pair of indices of the nodes of a resistor
  struct indexPair {
    /// Row of the first node
    std::size_t row1;
    /// Column of the first node
    std::size_t row2;
    /// Row of the second node
    std::size_t row2;
    /// Column of the second node
    std::size_t col2;
  };
}


#endif //PROYECTO2_RESISTORGRID_H
