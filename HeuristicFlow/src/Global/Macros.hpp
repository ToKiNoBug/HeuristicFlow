#ifndef HEU_MACROS_HPP
#define HEU_MACROS_HPP

#define HEU_RELOAD_MEMBERFUCTION_RUN \
  void run() { this->template __impl_run<typename std::decay<decltype(*this)>::type>(); }

#endif  //  HEU_MACROS_HPP