/*******************************************************************************
* Author    :  Angus Johnson                                                   *
* Date      :  12 October 2025                                                 *
* Website   :  https://www.angusj.com                                          *
* Copyright :  Angus Johnson 2010-2025                                         *
* Purpose   :  Core Clipper Library structures and functions                   *
* License   :  https://www.boost.org/LICENSE_1_0.txt                           *
*******************************************************************************/

#ifndef CLIPPER_CORE_H
#define CLIPPER_CORE_H

#include "clipper2/clipper.version.h"
#include <cstdint>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <cmath>

namespace Clipper2Lib
{

#if (defined(__cpp_exceptions) && __cpp_exceptions) || (defined(__EXCEPTIONS) && __EXCEPTIONS)

  class Clipper2Exception : public std::exception {
  public:
    explicit Clipper2Exception(const char* description) :
      m_descr(description) {}
    virtual const char* what() const noexcept override { return m_descr.c_str(); }
  private:
    std::string m_descr;
  };

  static const char* precision_error =
    "Precision exceeds the permitted range";
  static const char* range_error =
    "Values exceed permitted range";
  static const char* scale_error =
    "Invalid scale (either 0 or too large)";
  static const char* non_pair_error =
    "There must be 2 values for each coordinate";
  static const char* undefined_error =
    "There is an undefined error in Clipper2";
#endif

  // error codes (2^n)
  const int precision_error_i   = 1;  // non-fatal
  const int scale_error_i       = 2;  // non-fatal
  const int non_pair_error_i    = 4;  // non-fatal
  const int undefined_error_i   = 32; // fatal
  const int range_error_i       = 64;

#ifndef PI
  static const double PI = 3.141592653589793238;
#endif

#ifdef CLIPPER2_MAX_DECIMAL_PRECISION
  const int CLIPPER2_MAX_DEC_PRECISION = CLIPPER2_MAX_DECIMAL_PRECISION;
#else
  const int CLIPPER2_MAX_DEC_PRECISION = 8; // see Discussions #564
#endif

  static const int64_t MAX_COORD = INT64_MAX >> 2;
  static const int64_t MIN_COORD = -MAX_COORD;
  static const int64_t INVALID = INT64_MAX;
  const double max_coord = static_cast<double>(MAX_COORD);
  const double min_coord = static_cast<double>(MIN_COORD);

  static const double MAX_DBL = (std::numeric_limits<double>::max)();

  static void DoError([[maybe_unused]] int error_code)
  {
#if (defined(__cpp_exceptions) && __cpp_exceptions) || (defined(__EXCEPTIONS) && __EXCEPTIONS)
    switch (error_code)
    {
    case precision_error_i:
      throw Clipper2Exception(precision_error);
    case scale_error_i:
      throw Clipper2Exception(scale_error);
    case non_pair_error_i:
      throw Clipper2Exception(non_pair_error);
    case undefined_error_i:
      throw Clipper2Exception(undefined_error);
    case range_error_i:
      throw Clipper2Exception(range_error);
    // Should never happen, but adding this to stop a compiler warning
    default:
      throw Clipper2Exception("Unknown error");
    }
#else
    ++error_code; // only to stop compiler warning
#endif
  }

  // can we call std::round on T? (default false) (#824)
  template <typename T, typename = void>
  struct is_round_invocable : std::false_type {};

  template <typename T>
  struct is_round_invocable<T, std::void_t<decltype(std::round(std::declval<T>()))>> : std::true_type {};


  //By far the most widely used filling rules for polygons are EvenOdd
  //and NonZero, sometimes called Alternate and Winding respectively.
  //https://en.wikipedia.org/wiki/Nonzero-rule
  enum class FillRule { EvenOdd, NonZero, Positive, Negative };

#ifdef USINGZ
  using z_type = int64_t;
#endif

  // Point ------------------------------------------------------------------------

  template <typename T>
  struct Point {
    T x;
    T y;
#ifdef USINGZ
     z_type z;

    template <typename T2>
    inline void Init(const T2 x_ = 0, const T2 y_ = 0, const z_type z_ = 0)
    {
      if constexpr (std::is_integral_v<T> &&
        is_round_invocable<T2>::value && !std::is_integral_v<T2>)
      {
        x = static_cast<T>(std::round(x_));
        y = static_cast<T>(std::round(y_));
        z = z_;
      }
      else
      {
        x = static_cast<T>(x_);
        y = static_cast<T>(y_);
        z = z_;
      }
    }

    explicit Point() : x(0), y(0), z(0) {}

    template <typename T2>
    Point(const T2 x_, const T2 y_, const z_type z_ = 0)
    {
      Init(x_, y_);
      z = z_;
    }

    template <typename T2>
    explicit Point(const Point<T2>& p)
    {
      Init(p.x, p.y, p.z);
    }

    template <typename T2>
    explicit Point(const Point<T2>& p, z_type z_)
    {
      Init(p.x, p.y, z_);
    }

    Point operator * (const double scale) const
    {
      return Point(x * scale, y * scale, z);
    }

    void SetZ(const z_type z_value) { z = z_value; }

    friend std::ostream& operator<<(std::ostream& os, const Point& point)
    {
      os << point.x << "," << point.y << "," << point.z;
      return os;
    }

#else

    template <typename T2>
    inline void Init(const T2 x_ = 0, const T2 y_ = 0)
    {
      if constexpr (std::is_integral_v<T> &&
        is_round_invocable<T2>::value && !std::is_integral_v<T2>)
      {
        x = static_cast<T>(std::round(x_));
        y = static_cast<T>(std::round(y_));
      }
      else
      {
        x = static_cast<T>(x_);
        y = static_cast<T>(y_);
      }
    }

    explicit Point() : x(0), y(0) {}

    template <typename T2>
    Point(const T2 x_, const T2 y_) { Init(x_, y_); }

    template <typename T2>
    explicit Point(const Point<T2>& p) { Init(p.x, p.y); }

    Point operator * (const double scale) const
    {
      return Point(x * scale, y * scale);
    }

    friend std::ostream& operator<<(std::ostream& os, const Point& point)
    {
      os << point.x << "," << point.y;
      return os;
    }
#endif

    friend bool operator==(const Point& a, const Point& b)
    {
      return a.x == b.x && a.y == b.y;
    }

    friend bool operator!=(const Point& a, const Point& b)
    {
      return !(a == b);
    }

    inline Point<T> operator-() const
    {
      return Point<T>(-x, -y);
    }

    inline Point operator+(const Point& b) const
    {
      return Point(x + b.x, y + b.y);
    }

    inline Point operator-(const Point& b) const
    {
      return Point(x - b.x, y - b.y);
    }

    inline void Negate() { x = -x; y = -y; }

  };

  //nb: using 'using' here (instead of typedef) as they can be used in templates
  using Point64 = Point<int64_t>;
  using PointD = Point<double>;

  template <typename T>
  using Path = std::vector<Point<T>>;
  template <typename T>
  using Paths = std::vector<Path<T>>;

  template <typename T, typename T2=T>
  Path<T>& operator<<(Path<T>& poly, const Point<T2>& p)
  {
    poly.emplace_back(p);
    return poly;
  }

  template <typename T>
  Paths<T>& operator<<(Paths<T>& polys, const Path<T>& p)
  {
    polys.emplace_back(p);
    return polys;
  }

  using Path64 = Path<int64_t>;
  using PathD = Path<double>;
  using Paths64 = std::vector< Path64>;
  using PathsD = std::vector< PathD>;

  static const Point64 InvalidPoint64 = Point64(
    (std::numeric_limits<int64_t>::max)(),
    (std::numeric_limits<int64_t>::max)());
  static const PointD InvalidPointD = PointD(
    (std::numeric_limits<double>::max)(),
    (std::numeric_limits<double>::max)());

  template<typename T>
  static inline Point<T> MidPoint(const Point<T>& p1, const Point<T>& p2)
  {
    Point<T> result;
    result.x = (p1.x + p2.x) / 2;
    result.y = (p1.y + p2.y) / 2;
    return result;
  }

  // Rect ------------------------------------------------------------------------

  template <typename T>
  struct Rect;

  using Rect64 = Rect<int64_t>;
  using RectD = Rect<double>;

  template <typename T>
  struct Rect {
    T left;
    T top;
    T right;
    T bottom;

    Rect(T l, T t, T r, T b) :
      left(l),
      top(t),
      right(r),
      bottom(b) {}

    Rect(bool is_valid = true)
    {
      if (is_valid)
      {
        left = right = top = bottom = 0;
      }
      else
      {
        left = top = (std::numeric_limits<T>::max)();
        right = bottom = std::numeric_limits<T>::lowest();
      }
    }

    static Rect<T> InvalidRect()
    {
      return {
        (std::numeric_limits<T>::max)(),
        (std::numeric_limits<T>::max)(),
        std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest() };
    }

    bool IsValid() const { return left != (std::numeric_limits<T>::max)(); }

    T Width() const { return right - left; }
    T Height() const { return bottom - top; }
    void Width(T width) { right = left + width; }
    void Height(T height) { bottom = top + height; }

    Point<T> MidPoint() const
    {
      return Point<T>((left + right) / 2, (top + bottom) / 2);
    }

    Path<T> AsPath() const
    {
      Path<T> result;
      result.reserve(4);
      result.emplace_back(left, top);
      result.emplace_back(right, top);
      result.emplace_back(right, bottom);
      result.emplace_back(left, bottom);
      return result;
    }

    bool Contains(const Point<T>& pt) const
    {
      return pt.x > left && pt.x < right&& pt.y > top && pt.y < bottom;
    }

    bool Contains(const Rect<T>& rec) const
    {
      return rec.left >= left && rec.right <= right &&
        rec.top >= top && rec.bottom <= bottom;
    }

    void Scale(double scale) {
      left *= scale;
      top *= scale;
      right *= scale;
      bottom *= scale;
    }

    bool IsEmpty() const { return bottom <= top || right <= left; }

    bool Intersects(const Rect<T>& rec) const
    {
      return ((std::max)(left, rec.left) <= (std::min)(right, rec.right)) &&
        ((std::max)(top, rec.top) <= (std::min)(bottom, rec.bottom));
    }

    bool operator==(const Rect<T>& other) const {
      return left == other.left && right == other.right &&
        top == other.top && bottom == other.bottom;
    }

    Rect<T>& operator+=(const Rect<T>& other)
    {
      left = (std::min)(left, other.left);
      top = (std::min)(top, other.top);
      right = (std::max)(right, other.right);
      bottom = (std::max)(bottom, other.bottom);
      return *this;
    }

    Rect<T> operator+(const Rect<T>& other) const
    {
      Rect<T> result = *this;
      result += other;
      return result;
    }

    friend std::ostream& operator<<(std::ostream& os, const Rect<T>& rect) {
      os << "(" << rect.left << "," << rect.top << "," << rect.right << "," << rect.bottom << ") ";
      return os;
    }
  };

  template <typename T1, typename T2>
  inline Rect<T1> ScaleRect(const Rect<T2>& rect, double scale)
  {
    Rect<T1> result;

    if constexpr (std::is_integral_v<T1> &&
      is_round_invocable<T2>::value && !std::is_integral_v<T2>)
    {
      result.left = static_cast<T1>(std::round(rect.left * scale));
      result.top = static_cast<T1>(std::round(rect.top * scale));
      result.right = static_cast<T1>(std::round(rect.right * scale));
      result.bottom = static_cast<T1>(std::round(rect.bottom * scale));
    }
    else
    {
      result.left = static_cast<T1>(rect.left * scale);
      result.top = static_cast<T1>(rect.top * scale);
      result.right = static_cast<T1>(rect.right * scale);
      result.bottom = static_cast<T1>(rect.bottom * scale);
    }
    return result;
  }

  static const Rect64 InvalidRect64 = Rect64::InvalidRect();
  static const RectD InvalidRectD = RectD::InvalidRect();

  template <typename T>
  Rect<T> GetBounds(const Path<T>& path)
  {
    T xmin = (std::numeric_limits<T>::max)();
    T ymin = (std::numeric_limits<T>::max)();
    T xmax = std::numeric_limits<T>::lowest();
    T ymax = std::numeric_limits<T>::lowest();
    for (const auto& p : path)
    {
      if (p.x < xmin) xmin = p.x;
      if (p.x > xmax) xmax = p.x;
      if (p.y < ymin) ymin = p.y;
      if (p.y > ymax) ymax = p.y;
    }
    return Rect<T>(xmin, ymin, xmax, ymax);
  }

  template <typename T>
  Rect<T> GetBounds(const Paths<T>& paths)
  {
    T xmin = (std::numeric_limits<T>::max)();
    T ymin = (std::numeric_limits<T>::max)();
    T xmax = std::numeric_limits<T>::lowest();
    T ymax = std::numeric_limits<T>::lowest();
    for (const Path<T>& path : paths)
      for (const Point<T>& p : path)
      {
        if (p.x < xmin) xmin = p.x;
        if (p.x > xmax) xmax = p.x;
        if (p.y < ymin) ymin = p.y;
        if (p.y > ymax) ymax = p.y;
      }
    return Rect<T>(xmin, ymin, xmax, ymax);
  }

  template <typename T, typename T2>
  Rect<T> GetBounds(const Path<T2>& path)
  {
    T xmin = (std::numeric_limits<T>::max)();
    T ymin = (std::numeric_limits<T>::max)();
    T xmax = std::numeric_limits<T>::lowest();
    T ymax = std::numeric_limits<T>::lowest();
    for (const auto& p : path)
    {
      if (p.x < xmin) xmin = static_cast<T>(p.x);
      if (p.x > xmax) xmax = static_cast<T>(p.x);
      if (p.y < ymin) ymin = static_cast<T>(p.y);
      if (p.y > ymax) ymax = static_cast<T>(p.y);
    }
    return Rect<T>(xmin, ymin, xmax, ymax);
  }

  template <typename T, typename T2>
  Rect<T> GetBounds(const Paths<T2>& paths)
  {
    T xmin = (std::numeric_limits<T>::max)();
    T ymin = (std::numeric_limits<T>::max)();
    T xmax = std::numeric_limits<T>::lowest();
    T ymax = std::numeric_limits<T>::lowest();
    for (const Path<T2>& path : paths)
      for (const Point<T2>& p : path)
      {
        if (p.x < xmin) xmin = static_cast<T>(p.x);
        if (p.x > xmax) xmax = static_cast<T>(p.x);
        if (p.y < ymin) ymin = static_cast<T>(p.y);
        if (p.y > ymax) ymax = static_cast<T>(p.y);
      }
    return Rect<T>(xmin, ymin, xmax, ymax);
  }

  template <typename T>
  std::ostream& operator << (std::ostream& outstream, const Path<T>& path)
  {
    if (!path.empty())
    {
      auto pt = path.cbegin(), last = path.cend() - 1;
      while (pt != last)
        outstream << *pt++ << ", ";
      outstream << *last << std::endl;
    }
    return outstream;
  }

  template <typename T>
  std::ostream& operator << (std::ostream& outstream, const Paths<T>& paths)
  {
    for (auto p : paths)
      outstream << p;
    return outstream;
  }


  template <typename T1, typename T2>
  inline Path<T1> ScalePath(const Path<T2>& path,
    double scale_x, double scale_y, int& error_code)
  {
    Path<T1> result;
    if (scale_x == 0 || scale_y == 0)
    {
      error_code |= scale_error_i;
      DoError(scale_error_i);
      // if no exception, treat as non-fatal error
      if (scale_x == 0) scale_x = 1.0;
      if (scale_y == 0) scale_y = 1.0;
    }

    result.reserve(path.size());
#ifdef USINGZ
    std::transform(path.begin(), path.end(), back_inserter(result),
      [scale_x, scale_y](const auto& pt)
      { return Point<T1>(pt.x * scale_x, pt.y * scale_y, pt.z); });
#else
    std::transform(path.begin(), path.end(), back_inserter(result),
      [scale_x, scale_y](const auto& pt)
      { return Point<T1>(pt.x * scale_x, pt.y * scale_y); });
#endif
    return result;
  }

  template <typename T1, typename T2>
  inline Path<T1> ScalePath(const Path<T2>& path,
    double scale, int& error_code)
  {
    return ScalePath<T1, T2>(path, scale, scale, error_code);
  }

  template <typename T1, typename T2>
  inline Paths<T1> ScalePaths(const Paths<T2>& paths,
    double scale_x, double scale_y, int& error_code)
  {
    Paths<T1> result;

    if constexpr (std::is_integral_v<T1>)
    {
      RectD r = GetBounds<double, T2>(paths);
      if ((r.left * scale_x) < min_coord ||
        (r.right * scale_x) > max_coord ||
        (r.top * scale_y) < min_coord ||
        (r.bottom * scale_y) > max_coord)
      {
        error_code |= range_error_i;
        DoError(range_error_i);
        return result; // empty path
      }
    }

    result.reserve(paths.size());
    std::transform(paths.begin(), paths.end(), back_inserter(result),
      [=, &error_code](const auto& path)
      { return ScalePath<T1, T2>(path, scale_x, scale_y, error_code); });
    return result;
  }

  template <typename T1, typename T2>
  inline Paths<T1> ScalePaths(const Paths<T2>& paths,
    double scale, int& error_code)
  {
    return ScalePaths<T1, T2>(paths, scale, scale, error_code);
  }

  template <typename T1, typename T2>
  inline Path<T1> TransformPath(const Path<T2>& path)
  {
    Path<T1> result;
    result.reserve(path.size());
    std::transform(path.cbegin(), path.cend(), std::back_inserter(result),
      [](const Point<T2>& pt) {return Point<T1>(pt); });
    return result;
  }

  template <typename T1, typename T2>
  inline Paths<T1> TransformPaths(const Paths<T2>& paths)
  {
    Paths<T1> result;
    std::transform(paths.cbegin(), paths.cend(), std::back_inserter(result),
      [](const Path<T2>& path) {return TransformPath<T1, T2>(path); });
    return result;
  }

  template<typename T>
  inline double Sqr(T val)
  {
    return static_cast<double>(val) * static_cast<double>(val);
  }

  template<typename T>
  inline bool NearEqual(const Point<T>& p1,
    const Point<T>& p2, double max_dist_sqrd)
  {
    return Sqr(p1.x - p2.x) + Sqr(p1.y - p2.y) < max_dist_sqrd;
  }

  template<typename T>
  inline Path<T> StripNearEqual(const Path<T>& path,
    double max_dist_sqrd, bool is_closed_path)
  {
    if (path.size() == 0) return Path<T>();
    Path<T> result;
    result.reserve(path.size());
    typename Path<T>::const_iterator path_iter = path.cbegin();
    Point<T> first_pt = *path_iter++, last_pt = first_pt;
    result.emplace_back(first_pt);
    for (; path_iter != path.cend(); ++path_iter)
    {
      if (!NearEqual(*path_iter, last_pt, max_dist_sqrd))
      {
        last_pt = *path_iter;
        result.emplace_back(last_pt);
      }
    }
    if (!is_closed_path) return result;
    while (result.size() > 1 &&
      NearEqual(result.back(), first_pt, max_dist_sqrd)) result.pop_back();
    return result;
  }

  template<typename T>
  inline Paths<T> StripNearEqual(const Paths<T>& paths,
    double max_dist_sqrd, bool is_closed_path)
  {
    Paths<T> result;
    result.reserve(paths.size());
    for (typename Paths<T>::const_iterator paths_citer = paths.cbegin();
      paths_citer != paths.cend(); ++paths_citer)
    {
      result.emplace_back(std::move(StripNearEqual(*paths_citer, max_dist_sqrd, is_closed_path)));
    }
    return result;
  }

  template<typename T>
  inline void StripDuplicates( Path<T>& path, bool is_closed_path)
  {
    //https://stackoverflow.com/questions/1041620/whats-the-most-efficient-way-to-erase-duplicates-and-sort-a-vector#:~:text=Let%27s%20compare%20three%20approaches%3A
    path.erase(std::unique(path.begin(), path.end()), path.end());
    if (is_closed_path)
      while (path.size() > 1 && path.back() == path.front()) path.pop_back();
  }

  template<typename T>
  inline void StripDuplicates( Paths<T>& paths, bool is_closed_path)
  {
    for (typename Paths<T>::iterator paths_citer = paths.begin();
      paths_citer != paths.end(); ++paths_citer)
    {
      StripDuplicates(*paths_citer, is_closed_path);
    }
  }

  // Miscellaneous ------------------------------------------------------------

  inline void CheckPrecisionRange(int& precision, int& error_code)
  {
    if (precision >= -CLIPPER2_MAX_DEC_PRECISION &&
      precision <= CLIPPER2_MAX_DEC_PRECISION) return;
    error_code |= precision_error_i; // non-fatal error
    DoError(precision_error_i);      // does nothing when exceptions are disabled
    precision = precision > 0 ? CLIPPER2_MAX_DEC_PRECISION : -CLIPPER2_MAX_DEC_PRECISION;
  }

  inline void CheckPrecisionRange(int& precision)
  {
    int error_code = 0;
    CheckPrecisionRange(precision, error_code);
  }

  inline int TriSign(int64_t x) // returns 0, 1 or -1
  {
    return (x > 0) - (x < 0); 
  }

  struct UInt128Struct
  {
    const uint64_t lo = 0;
    const uint64_t hi = 0;

    bool operator==(const UInt128Struct& other) const
    {
      return lo == other.lo && hi == other.hi;
    }

  };

  inline UInt128Struct MultiplyUInt64(uint64_t a, uint64_t b) // #834, #835
  {
    // note to self - lamba expressions follow
    const auto lo = [](uint64_t x) { return x & 0xFFFFFFFF; };
    const auto hi = [](uint64_t x) { return x >> 32; };

    const uint64_t x1 = lo(a) * lo(b);
    const uint64_t x2 = hi(a) * lo(b) + hi(x1);
    const uint64_t x3 = lo(a) * hi(b) + lo(x2);
    return { uint64_t(lo(x3) << 32 | lo(x1)), uint64_t(hi(a) * hi(b) + hi(x2) + hi(x3)) };
  }

  // returns true if (and only if) a * b == c * d
  inline bool ProductsAreEqual(int64_t a, int64_t b, int64_t c, int64_t d)
  {
#if (defined(__clang__) || defined(__GNUC__)) && UINTPTR_MAX >= UINT64_MAX
    const auto ab = static_cast<__int128_t>(a) * static_cast<__int128_t>(b);
    const auto cd = static_cast<__int128_t>(c) * static_cast<__int128_t>(d);
    return ab == cd;
#else
    // nb: unsigned values needed for calculating overflow carry
    const auto abs_a = static_cast<uint64_t>(std::abs(a));
    const auto abs_b = static_cast<uint64_t>(std::abs(b));
    const auto abs_c = static_cast<uint64_t>(std::abs(c));
    const auto abs_d = static_cast<uint64_t>(std::abs(d));

    const auto ab = MultiplyUInt64(abs_a, abs_b);
    const auto cd = MultiplyUInt64(abs_c, abs_d);

    // nb: it's important to differentiate 0 values here from other values
    const auto sign_ab = TriSign(a) * TriSign(b);
    const auto sign_cd = TriSign(c) * TriSign(d);

    return ab == cd && sign_ab == sign_cd;
#endif
  }

  template <typename T>
  inline int CrossProductSign(const Point<T>& pt1, const Point<T>& pt2, const Point<T>& pt3)
  {
    const auto a = pt2.x - pt1.x;
    const auto b = pt3.y - pt2.y;
    const auto c = pt2.y - pt1.y;
    const auto d = pt3.x - pt2.x;

#if (defined(__clang__) || defined(__GNUC__)) && UINTPTR_MAX >= UINT64_MAX
    const auto ab = static_cast<__int128_t>(a) * static_cast<__int128_t>(b);
    const auto cd = static_cast<__int128_t>(c) * static_cast<__int128_t>(d);
    if (ab > cd) return 1;
    else if (ab < cd) return -1;
    else return 0;
#else
    const auto ab = MultiplyUInt64(std::abs(a), std::abs(b));
    const auto cd = MultiplyUInt64(std::abs(c), std::abs(d));

    const auto sign_ab = TriSign(a) * TriSign(b);
    const auto sign_cd = TriSign(c) * TriSign(d);

    if (sign_ab == sign_cd)
    {
      int result;
      if (ab.hi == cd.hi)
      {
        if (ab.lo == cd.lo) return 0;
        result = (ab.lo > cd.lo) ? 1 : -1;
      }
      else result = (ab.hi > cd.hi) ? 1 : -1;
      return (sign_ab > 0) ? result : -result;
    }
    return (sign_ab > sign_cd) ? 1 : -1;
#endif
  }

  template <typename T>
  inline bool IsCollinear(const Point<T>& pt1,
    const Point<T>& sharedPt, const Point<T>& pt2) // #777
  {
    const auto a = sharedPt.x - pt1.x;
    const auto b = pt2.y - sharedPt.y;
    const auto c = sharedPt.y - pt1.y;
    const auto d = pt2.x - sharedPt.x;
    // When checking for collinearity with very large coordinate values
    // then ProductsAreEqual is more accurate than using CrossProduct.
    return ProductsAreEqual(a, b, c, d);
  }


  template <typename T>
  inline double CrossProduct(const Point<T>& pt1, const Point<T>& pt2, const Point<T>& pt3) {
    return (static_cast<double>(pt2.x - pt1.x) * static_cast<double>(pt3.y -
      pt2.y) - static_cast<double>(pt2.y - pt1.y) * static_cast<double>(pt3.x - pt2.x));
  }

  template <typename T>
  inline double CrossProduct(const Point<T>& vec1, const Point<T>& vec2)
  {
    return static_cast<double>(vec1.y * vec2.x) - static_cast<double>(vec2.y * vec1.x);
  }

  template <typename T>
  inline double DotProduct(const Point<T>& pt1, const Point<T>& pt2, const Point<T>& pt3) {
    return (static_cast<double>(pt2.x - pt1.x) * static_cast<double>(pt3.x - pt2.x) +
      static_cast<double>(pt2.y - pt1.y) * static_cast<double>(pt3.y - pt2.y));
  }

  template <typename T>
  inline double DotProduct(const Point<T>& vec1, const Point<T>& vec2)
  {
    return static_cast<double>(vec1.x * vec2.x) + static_cast<double>(vec1.y * vec2.y);
  }

  template <typename T>
  inline double DistanceSqr(const Point<T> pt1, const Point<T> pt2)
  {
    return Sqr(pt1.x - pt2.x) + Sqr(pt1.y - pt2.y);
  }

  template <typename T>
  inline double PerpendicDistFromLineSqrd(const Point<T>& pt,
    const Point<T>& line1, const Point<T>& line2)
  {
    //perpendicular distance of point (x³,y³) = (Ax³ + By³ + C)/Sqrt(A² + B²)
    //see https://en.wikipedia.org/wiki/Perpendicular_distance
    double a = static_cast<double>(pt.x - line1.x);
    double b = static_cast<double>(pt.y - line1.y);
    double c = static_cast<double>(line2.x - line1.x);
    double d = static_cast<double>(line2.y - line1.y);
    if (c == 0 && d == 0) return 0;
    return Sqr(a * d - c * b) / (c * c + d * d);
  }

  template <typename T>
  inline double Area(const Path<T>& path)
  {
    size_t cnt = path.size();
    if (cnt < 3) return 0.0;
    double a = 0.0;
    typename Path<T>::const_iterator it1, it2 = path.cend() - 1, stop = it2;
    if (!(cnt & 1)) ++stop;
    for (it1 = path.cbegin(); it1 != stop;)
    {
      a += static_cast<double>(it2->y + it1->y) * (it2->x - it1->x);
      it2 = it1 + 1;
      a += static_cast<double>(it1->y + it2->y) * (it1->x - it2->x);
      it1 += 2;
    }
    if (cnt & 1)
      a += static_cast<double>(it2->y + it1->y) * (it2->x - it1->x);
    return (a * 0.5);
  }

  template <typename T>
  inline double Area(const Paths<T>& paths)
  {
    double a = 0.0;
    for (typename Paths<T>::const_iterator paths_iter = paths.cbegin();
      paths_iter != paths.cend(); ++paths_iter)
    {
      a += Area<T>(*paths_iter);
    }
    return a;
  }

  template <typename T>
  inline bool IsPositive(const Path<T>& poly)
  {
    // A curve has positive orientation [and area] if a region 'R'
    // is on the left when traveling around the outside of 'R'.
    //https://mathworld.wolfram.com/CurveOrientation.html
    //nb: This statement is premised on using Cartesian coordinates
    return Area<T>(poly) >= 0;
  }

  // --------------------------------------------------------------------------
  // Path statistics, self-intersection detection and area conservation checks.
  // These are read-only helpers: no clipping is performed.
  // --------------------------------------------------------------------------

  struct PathStats
  {
    size_t path_count = 0;    // number of paths included in the statistics
    size_t vertex_count = 0;  // total number of vertices
    double area = 0.0;        // signed area (sum over all paths)
    double abs_area = 0.0;    // sum of the absolute area of each path
    double perimeter = 0.0;   // total perimeter of all paths
    bool is_positive = true;  // signed area >= 0
  };

  struct AreaCheckResult
  {
    double expected_area = 0.0;
    double actual_area = 0.0;
    double abs_error = 0.0;
    double rel_error = 0.0;
    double tolerance = 0.0;
    bool passed = false;
  };

  template <typename T>
  inline double Perimeter(const Path<T>& path, bool is_closed_path = true)
  {
    const size_t cnt = path.size();
    if (cnt < 2) return 0.0;
    double result = 0.0;
    for (size_t i = 1; i < cnt; ++i)
      result += std::sqrt(Sqr(path[i].x - path[i - 1].x) +
        Sqr(path[i].y - path[i - 1].y));
    if (is_closed_path)
      result += std::sqrt(Sqr(path[0].x - path[cnt - 1].x) +
        Sqr(path[0].y - path[cnt - 1].y));
    return result;
  }

  template <typename T>
  inline double Perimeter(const Paths<T>& paths, bool is_closed_path = true)
  {
    double result = 0.0;
    for (const Path<T>& path : paths)
      result += Perimeter(path, is_closed_path);
    return result;
  }

  template <typename T>
  inline size_t VertexCount(const Path<T>& path)
  {
    return path.size();
  }

  template <typename T>
  inline size_t VertexCount(const Paths<T>& paths)
  {
    size_t result = 0;
    for (const Path<T>& path : paths) result += path.size();
    return result;
  }

  template <typename T>
  inline PathStats GetPathStats(const Path<T>& path,
    bool is_closed_path = true)
  {
    PathStats result;
    result.path_count = 1;
    result.vertex_count = path.size();
    result.area = Area(path);
    result.abs_area = std::fabs(result.area);
    result.perimeter = Perimeter(path, is_closed_path);
    result.is_positive = result.area >= 0.0;
    return result;
  }

  template <typename T>
  inline PathStats GetPathStats(const Paths<T>& paths,
    bool is_closed_path = true)
  {
    PathStats result;
    result.path_count = paths.size();
    result.area = Area(paths);
    for (const Path<T>& path : paths)
    {
      result.vertex_count += path.size();
      result.abs_area += std::fabs(Area(path));
      result.perimeter += Perimeter(path, is_closed_path);
    }
    result.is_positive = result.area >= 0.0;
    return result;
  }

  namespace detail
  {
    // Unlike CrossProductSign (which truncates non-integral arguments to
    // 128 bit integers), OrientSign also handles floating point points.
    template <typename T>
    inline int OrientSign(const Point<T>& pt1,
      const Point<T>& pt2, const Point<T>& pt3)
    {
      if constexpr (std::is_integral_v<T>)
        return CrossProductSign(pt1, pt2, pt3);
      else
      {
        const T ab = (pt2.x - pt1.x) * (pt3.y - pt2.y);
        const T cd = (pt2.y - pt1.y) * (pt3.x - pt2.x);
        return (ab > cd) ? 1 : ((ab < cd) ? -1 : 0);
      }
    }

    template <typename T>
    inline bool PointInBoundingBox(const Point<T>& pt,
      const Point<T>& pt1, const Point<T>& pt2)
    {
      return pt.x >= (std::min)(pt1.x, pt2.x) &&
        pt.x <= (std::max)(pt1.x, pt2.x) &&
        pt.y >= (std::min)(pt1.y, pt2.y) &&
        pt.y <= (std::max)(pt1.y, pt2.y);
    }

    // Returns true when the closed segments a-b and c-d share at least
    // one point (crossings, endpoint touches and collinear overlaps).
    template <typename T>
    inline bool SegmentsTouch(const Point<T>& a, const Point<T>& b,
      const Point<T>& c, const Point<T>& d)
    {
      const int o1 = OrientSign(a, b, c);
      const int o2 = OrientSign(a, b, d);
      const int o3 = OrientSign(c, d, a);
      const int o4 = OrientSign(c, d, b);
      if (o1 != o2 && o3 != o4) return true;
      if (o1 == 0 && PointInBoundingBox(c, a, b)) return true;
      if (o2 == 0 && PointInBoundingBox(d, a, b)) return true;
      if (o3 == 0 && PointInBoundingBox(a, c, d)) return true;
      if (o4 == 0 && PointInBoundingBox(b, c, d)) return true;
      return false;
    }

    // Consecutive edges (pt1, shared) and (shared, pt2) should only meet
    // at 'shared'. They also self-overlap when they collinearly double
    // back over each other (eg a spike: A -> B -> A).
    template <typename T>
    inline bool AdjacentEdgesOverlap(const Point<T>& pt1,
      const Point<T>& shared, const Point<T>& pt2)
    {
      if (OrientSign(pt1, shared, pt2) != 0) return false;
      const double dx1 = static_cast<double>(pt1.x) - shared.x;
      const double dy1 = static_cast<double>(pt1.y) - shared.y;
      const double dx2 = static_cast<double>(pt2.x) - shared.x;
      const double dy2 = static_cast<double>(pt2.y) - shared.y;
      return (dx1 * dx2 + dy1 * dy2) > 0.0;
    }
  } // namespace detail

  // Counts pairs of edges in a closed path that overlap or cross. Edges that
  // merely share a vertex with their immediate neighbour are not counted,
  // but adjacent edges that retrace each other are. Zero length edges are
  // ignored (so paths that repeat their first point at the end are fine).
  // Intersections between different paths are not considered.
  template <typename T>
  inline size_t CountSelfIntersections(const Path<T>& path)
  {
    // Drop repeated consecutive vertices, and a repeated final vertex, so
    // that "adjacent" edges really are edge i and edge i+1 (cyclically).
    Path<T> pts;
    pts.reserve(path.size());
    for (const Point<T>& pt : path)
      if (pts.empty() || pts.back() != pt) pts.push_back(pt);
    if (pts.size() > 1 && pts.front() == pts.back()) pts.pop_back();

    const size_t n = pts.size();
    if (n < 3) return 0;

    struct Edge
    {
      size_t i0, i1;
      T minx, maxx, miny, maxy;
    };

    std::vector<Edge> edges;
    edges.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
      const size_t j = (i + 1 == n) ? 0 : i + 1;
      const Point<T>& a = pts[i];
      const Point<T>& b = pts[j];
      edges.push_back(Edge{ i, j,
        (std::min)(a.x, b.x), (std::max)(a.x, b.x),
        (std::min)(a.y, b.y), (std::max)(a.y, b.y) });
    }

    std::sort(edges.begin(), edges.end(),
      [](const Edge& a, const Edge& b)
      {
        if (a.minx != b.minx) return a.minx < b.minx;
        return a.maxx < b.maxx;
      });

    // Sweep a vertical line left to right, keeping an "active" set of edges
    // whose x extent is still open. With the edge list sorted by its left
    // x coordinate, each edge is only compared against still-active edges
    // whose left x is <= its own, so each pair is tested exactly once.
    size_t count = 0;
    std::vector<size_t> active;    // indices into edges
    std::vector<size_t> survivors;

    const auto test_pair = [&](size_t prev_idx, size_t cur_idx)
    {
      const Edge& prev = edges[prev_idx];
      const Edge& cur = edges[cur_idx];
      if (prev.maxy < cur.miny || prev.miny > cur.maxy) return;

      size_t shared = n; // n is not a valid vertex index
      if (prev.i0 == cur.i0 || prev.i0 == cur.i1) shared = prev.i0;
      else if (prev.i1 == cur.i0 || prev.i1 == cur.i1) shared = prev.i1;

      if (shared == n)
      {
        if (detail::SegmentsTouch(pts[cur.i0], pts[cur.i1],
          pts[prev.i0], pts[prev.i1])) ++count;
      }
      else
      {
        const size_t cur_other = (cur.i0 == shared) ? cur.i1 : cur.i0;
        const size_t prev_other = (prev.i0 == shared) ? prev.i1 : prev.i0;
        if (detail::AdjacentEdgesOverlap(pts[prev_other],
          pts[shared], pts[cur_other])) ++count;
      }
    };

    for (size_t e = 0; e < edges.size(); ++e)
    {
      // retire edges that end before this one begins
      const T x = edges[e].minx;
      survivors.clear();
      for (const size_t idx : active)
        if (edges[idx].maxx >= x) survivors.push_back(idx);
      active.swap(survivors);

      for (const size_t idx : active)
        test_pair(idx, e);

      active.push_back(e);
    }
    return count;
  }

  template <typename T>
  inline size_t CountSelfIntersections(const Paths<T>& paths)
  {
    size_t count = 0;
    for (const Path<T>& path : paths)
      count += CountSelfIntersections(path);
    return count;
  }

  template <typename T>
  inline bool IsSelfIntersecting(const Path<T>& path)
  {
    return CountSelfIntersections(path) > 0;
  }

  template <typename T>
  inline bool IsSelfIntersecting(const Paths<T>& paths)
  {
    return CountSelfIntersections(paths) > 0;
  }

  // Compares the (signed) area of 'output' against 'expected_area'.
  // 'tolerance' is an absolute area, so it's passed if
  // abs(actual - expected) <= tolerance.
  template <typename T>
  inline AreaCheckResult CheckAreaConservation(double expected_area,
    const Paths<T>& output, double tolerance)
  {
    AreaCheckResult result;
    result.expected_area = expected_area;
    result.actual_area = Area(output);
    result.tolerance = tolerance;
    result.abs_error =
      std::fabs(result.actual_area - result.expected_area);
    const double denom = std::fabs(result.expected_area);
    result.rel_error =
      (denom > 0.0) ? result.abs_error / denom : result.abs_error;
    result.passed = result.abs_error <= tolerance;
    return result;
  }

  // Compares the area of 'output' against the area of 'input'.
  template <typename T>
  inline AreaCheckResult CheckAreaConservation(const Paths<T>& input,
    const Paths<T>& output, double tolerance)
  {
    return CheckAreaConservation(Area(input), output, tolerance);
  }

  // GetLineIntersectPt - a 'true' result is non-parallel. The 'ip' will also
  // be constrained to seg1. However, it's possible that 'ip' won't be inside
  // seg2, even when 'ip' hasn't been constrained (ie 'ip' is inside seg1).

#if defined(CLIPPER2_HI_PRECISION) && CLIPPER2_HI_PRECISION
  // caution: this will compromise performance
  // https://github.com/AngusJohnson/Clipper2/issues/317#issuecomment-1314023253
  // See also CPP/BenchMark/GetIntersectPtBenchmark.cpp
  #define CC_MIN(x,y) ((x)>(y)?(y):(x))
  #define CC_MAX(x,y) ((x)<(y)?(y):(x))
  template<typename T>
  inline bool GetLineIntersectPt(const Point<T>& ln1a, const Point<T>& ln1b,
    const Point<T>& ln2a, const Point<T>& ln2b, Point<T>& ip)
  {
    double ln1dy = static_cast<double>(ln1b.y - ln1a.y);
    double ln1dx = static_cast<double>(ln1a.x - ln1b.x);
    double ln2dy = static_cast<double>(ln2b.y - ln2a.y);
    double ln2dx = static_cast<double>(ln2a.x - ln2b.x);
    double det = (ln2dy * ln1dx) - (ln1dy * ln2dx);
    if (det == 0.0) return false;
    T bb0minx = CC_MIN(ln1a.x, ln1b.x);
    T bb0miny = CC_MIN(ln1a.y, ln1b.y);
    T bb0maxx = CC_MAX(ln1a.x, ln1b.x);
    T bb0maxy = CC_MAX(ln1a.y, ln1b.y);
    T bb1minx = CC_MIN(ln2a.x, ln2b.x);
    T bb1miny = CC_MIN(ln2a.y, ln2b.y);
    T bb1maxx = CC_MAX(ln2a.x, ln2b.x);
    T bb1maxy = CC_MAX(ln2a.y, ln2b.y);

    if constexpr (std::is_integral_v<T>)
    {
      int64_t originx = (CC_MIN(bb0maxx, bb1maxx) + CC_MAX(bb0minx, bb1minx)) >> 1;
      int64_t originy = (CC_MIN(bb0maxy, bb1maxy) + CC_MAX(bb0miny, bb1miny)) >> 1;
      double ln0c = (ln1dy * static_cast<double>(ln1a.x - originx)) +
        (ln1dx * static_cast<double>(ln1a.y - originy));
      double ln1c = (ln2dy * static_cast<double>(ln2a.x - originx)) +
        (ln2dx * static_cast<double>(ln2a.y - originy));
      double hitx = ((ln1dx * ln1c) - (ln2dx * ln0c)) / det;
      double hity = ((ln2dy * ln0c) - (ln1dy * ln1c)) / det;

      ip.x = originx + (T)nearbyint(hitx);
      ip.y = originy + (T)nearbyint(hity);
    }
    else
    {
      double originx = (CC_MIN(bb0maxx, bb1maxx) + CC_MAX(bb0minx, bb1minx)) / 2.0;
      double originy = (CC_MIN(bb0maxy, bb1maxy) + CC_MAX(bb0miny, bb1miny)) / 2.0;
      double ln0c = (ln1dy * static_cast<double>(ln1a.x - originx)) +
        (ln1dx * static_cast<double>(ln1a.y - originy));
      double ln1c = (ln2dy * static_cast<double>(ln2a.x - originx)) +
        (ln2dx * static_cast<double>(ln2a.y - originy));
      double hitx = ((ln1dx * ln1c) - (ln2dx * ln0c)) / det;
      double hity = ((ln2dy * ln0c) - (ln1dy * ln1c)) / det;

      ip.x = originx + static_cast<T>(hitx);
      ip.y = originy + static_cast<T>(hity);
    }
#ifdef USINGZ
    ip.z = 0;
#endif
    return true;
}
#else
  template<typename T>
  inline bool GetLineIntersectPt(const Point<T>& ln1a, const Point<T>& ln1b,
    const Point<T>& ln2a, const Point<T>& ln2b, Point<T>& ip)
  {
    // https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    double dx1 = static_cast<double>(ln1b.x - ln1a.x);
    double dy1 = static_cast<double>(ln1b.y - ln1a.y);
    double dx2 = static_cast<double>(ln2b.x - ln2a.x);
    double dy2 = static_cast<double>(ln2b.y - ln2a.y);

    double det = dy1 * dx2 - dy2 * dx1;
    if (det == 0.0) return false;
    double t = ((ln1a.x - ln2a.x) * dy2 - (ln1a.y - ln2a.y) * dx2) / det;
    if (t <= 0.0) ip = ln1a;
    else if (t >= 1.0) ip = ln1b;
    else
    {
      ip.x = static_cast<T>(ln1a.x + t * dx1);
      ip.y = static_cast<T>(ln1a.y + t * dy1);
#ifdef USINGZ
      ip.z = 0;
#endif
    }
    return true;
  }
#endif

  template<typename T>
  inline Point<T> TranslatePoint(const Point<T>& pt, double dx, double dy)
  {
#ifdef USINGZ
    return Point<T>(pt.x + dx, pt.y + dy, pt.z);
#else
    return Point<T>(pt.x + dx, pt.y + dy);
#endif
  }


  template<typename T>
  inline Point<T> ReflectPoint(const Point<T>& pt, const Point<T>& pivot)
  {
#ifdef USINGZ
    return Point<T>(pivot.x + (pivot.x - pt.x), pivot.y + (pivot.y - pt.y), pt.z);
#else
    return Point<T>(pivot.x + (pivot.x - pt.x), pivot.y + (pivot.y - pt.y));
#endif
  }

  template<typename T>
  inline int GetSign(const T& val) 
  { 
    if (!val) return 0; 
    return (val > 0) ? 1 : -1;
  }

  inline bool SegmentsIntersect(const Point64& seg1a, const Point64& seg1b,
    const Point64& seg2a, const Point64& seg2b, bool inclusive = false)
  {
    double dy1 = static_cast<double>(seg1b.y - seg1a.y);
    double dx1 = static_cast<double>(seg1b.x - seg1a.x);
    double dy2 = static_cast<double>(seg2b.y - seg2a.y);
    double dx2 = static_cast<double>(seg2b.x - seg2a.x);
    double cp = dy1 * dx2 - dy2 * dx1;
    if (cp == 0) return false; // ie parallel segments

    if (inclusive)
    {
      //result **includes** segments that touch at an end point
      double t = ((seg1a.x - seg2a.x) * dy2 - (seg1a.y - seg2a.y) * dx2);
      if (t == 0) return true;
      if (t > 0)
      {
        if (cp < 0 || t > cp) return false;
      }
      else if (cp > 0 || t < cp) return false; // false when t more neg. than cp

      t = ((seg1a.x - seg2a.x) * dy1 - (seg1a.y - seg2a.y) * dx1);
      if (t == 0) return true;
      if (t > 0)  return (cp > 0 && t <= cp);
      else return (cp < 0 && t >= cp);        // true when t less neg. than cp
    }
    else 
    {
      //result **excludes** segments that touch at an end point
      double t = ((seg1a.x - seg2a.x) * dy2 - (seg1a.y - seg2a.y) * dx2);
      if (t == 0) return false;
      if (t > 0)
      {
        if (cp < 0 || t >= cp) return false;
      }
      else if (cp > 0 || t <= cp ) return false; // false when t more neg. than cp

      t = ((seg1a.x - seg2a.x) * dy1 - (seg1a.y - seg2a.y) * dx1);
      if (t == 0) return false;
      if (t > 0)  return (cp > 0 && t < cp);
      else return (cp < 0 && t > cp); // true when t less neg. than cp
    }
  }

  template<typename T>
  inline Point<T> GetClosestPointOnSegment(const Point<T>& offPt,
    const Point<T>& seg1, const Point<T>& seg2)
  {
    if (seg1.x == seg2.x && seg1.y == seg2.y) return seg1;
    double dx = static_cast<double>(seg2.x - seg1.x);
    double dy = static_cast<double>(seg2.y - seg1.y);
    double q =
      (static_cast<double>(offPt.x - seg1.x) * dx +
        static_cast<double>(offPt.y - seg1.y) * dy) /
      (Sqr(dx) + Sqr(dy));
    if (q < 0) q = 0; else if (q > 1) q = 1;
    if constexpr (std::is_integral_v<T>)
      return Point<T>(
        seg1.x + static_cast<T>(nearbyint(q * dx)),
        seg1.y + static_cast<T>(nearbyint(q * dy)));
    else
      return Point<T>(
        seg1.x + static_cast<T>(q * dx),
        seg1.y + static_cast<T>(q * dy));
  }

  enum class PointInPolygonResult { IsOn, IsInside, IsOutside };

  template <typename T>
  inline PointInPolygonResult PointInPolygon(const Point<T>& pt, const Path<T>& polygon)
  {
    if (polygon.size() < 3)
      return PointInPolygonResult::IsOutside;

    int val = 0;
    typename Path<T>::const_iterator cbegin = polygon.cbegin(), first = cbegin, curr, prev;
    typename Path<T>::const_iterator cend = polygon.cend();

    while (first != cend && first->y == pt.y) ++first;
    if (first == cend) // not a proper polygon
      return PointInPolygonResult::IsOutside;

    bool is_above = first->y < pt.y, starting_above = is_above;
    curr = first +1;
    while (true)
    {
      if (curr == cend)
      {
        if (cend == first || first == cbegin) break;
        cend = first;
        curr = cbegin;
      }

      if (is_above)
      {
        while (curr != cend && curr->y < pt.y) ++curr;
        if (curr == cend) continue;
      }
      else
      {
        while (curr != cend && curr->y > pt.y) ++curr;
        if (curr == cend) continue;
      }

      if (curr == cbegin)
        prev = polygon.cend() - 1; //nb: NOT cend (since might equal first)
      else
        prev = curr - 1;

      if (curr->y == pt.y)
      {
        if (curr->x == pt.x ||
          (curr->y == prev->y &&
            ((pt.x < prev->x) != (pt.x < curr->x))))
              return PointInPolygonResult::IsOn;
        ++curr;
        if (curr == first) break;
        continue;
      }

      if (pt.x < curr->x && pt.x < prev->x)
      {
        // we're only interested in edges crossing on the left
      }
      else if (pt.x > prev->x && pt.x > curr->x)
        val = 1 - val; // toggle val
      else
      {
        int d = CrossProductSign(*prev, *curr, pt);
        if (d == 0) return PointInPolygonResult::IsOn;
        if ((d < 0) == is_above) val = 1 - val;
      }
      is_above = !is_above;
      ++curr;
    }

    if (is_above != starting_above)
    {
      cend = polygon.cend();
      if (curr == cend) curr = cbegin;
      if (curr == cbegin) prev = cend - 1;
      else prev = curr - 1;
      int d = CrossProductSign(*prev, *curr, pt);
      if (d == 0) return PointInPolygonResult::IsOn;
      if ((d < 0) == is_above) val = 1 - val;
    }

    return (val == 0) ?
      PointInPolygonResult::IsOutside :
      PointInPolygonResult::IsInside;
  }

}  // namespace

#endif  // CLIPPER_CORE_H
