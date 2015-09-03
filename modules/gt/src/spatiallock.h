/*
    Copyright 2013 KU Leuven, Dept. Electrical Engineering, Medical Image Computing
    Herestraat 49 box 7003, 3000 Leuven, Belgium
    
    Written by Daan Christiaens, 17/03/14.
    
    This file is part of the Global Tractography module for MRtrix.
    
    MRtrix is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    MRtrix is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with MRtrix.  If not, see <http://www.gnu.org/licenses/>.
    
*/

#ifndef __gt_spatiallock_h__
#define __gt_spatiallock_h__

#include "point.h"
#include <mutex>

#include <set>


namespace MR {
  namespace DWI {
    namespace Tractography {
      namespace GT {

        /**
         * @brief SpatialLock manages a mutex lock on n positions in 3D space.
         */
        template <typename T = float >
        class SpatialLock
        {
        public:
          typedef T value_type;
          
          SpatialLock() : _tx(0), _ty(0), _tz(0) { }
          SpatialLock(const value_type t) : _tx(t), _ty(t), _tz(t) { }
          SpatialLock(const value_type tx, const value_type ty, const value_type tz) : _tx(tx), _ty(ty), _tz(tz) { }
          
          ~SpatialLock() {
            lockcentres.clear();
          }
          
          void setThreshold(const value_type t) {
            _tx = _ty = _tz = t;
          }
          
          void setThreshold(const value_type tx, const value_type ty, float tz) {
            _tx = tx;
            _ty = ty;
            _tz = tz;
          }
          
          bool lockIfNotLocked(const Point<T>& pos) {
            std::lock_guard<std::mutex> lock (mutex);
            Point<value_type> d;
            for (typename std::set<Point<value_type> >::iterator it = lockcentres.begin(); it != lockcentres.end(); ++it) {
              d = *it - pos;
              if ((std::abs(d[0]) < _tx) && (std::abs(d[1]) < _ty) && (std::abs(d[2]) < _tz))
                return false;
            }
            lockcentres.insert(pos);
            return true;
          }
          
          void unlock(const Point<T>& pos) {
            std::lock_guard<std::mutex> lock (mutex);
            lockcentres.erase(pos);
          }
          
        protected:
          std::mutex mutex;
          std::set<Point<value_type> > lockcentres;
          value_type _tx, _ty, _tz;
          
        };

      }
    }
  }
}

#endif // __gt_spatiallock_h__
