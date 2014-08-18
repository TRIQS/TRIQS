/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2014 by O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#pragma once
#include <vector>
#include <functional>

namespace triqs {
namespace utility {

 /**
  * A simple buffer for a generator.
  * Given a function, it provides a buffer of this function
  * Advantage :
  *  - do not pay the indirection cost at each call, but once every size call.
  *  - erase the function type
  * It is a semi-regular type.
  */
 template <typename R> struct buffered_function {

  /// Default constructor : no function bufferized. () will throw in this state
  buffered_function() = default;

  /** Constructor
   *
   * @tparam Function : type of the function to bufferize
   * @param f : function to bufferize
   * @param size : size of the buffer [optional]
   */
  template <typename Function> buffered_function(Function f, size_t size = 1000) : buffer(size) {
   refill = [f](buffered_function *bf) mutable { // without the mutable, the () of the lambda object is const, hence f
    for (auto &x : bf->buffer) x = f();
    bf->index = 0;
   };
   refill(this); // first filling of the buffer
  }

  /// Returns the next element. Refills the buffer if necessary.
  R operator()() {
   if (index > buffer.size() - 1) refill(this);
   return buffer[index++];
  }

  /// Returns the future next element, without increasing the index. Refills the buffer if necessary.
  R preview() {
   if (index > buffer.size() - 1) refill(this);
   return buffer[index];
  }

  private:
  size_t index;
  std::vector<R> buffer;
  std::function<void(buffered_function *)> refill; // this refills the buffer and reset index of a buffered_function.
  // NB : cannot capture this in refill because we want the object to be copyable and movable
 };
}
}
