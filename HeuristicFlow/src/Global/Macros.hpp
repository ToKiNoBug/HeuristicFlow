#ifndef HEU_MACROS_HPP
#define HEU_MACROS_HPP

#define HEU_RELOAD_MEMBERFUCTION_RUN                                         \
  inline void run() noexcept {                                               \
    this->template __impl_run<typename std::decay<decltype(*this)>::type>(); \
  }

#define HEU_DISPLINE \
  ::std::cout << "File : " << __FILE__ << " , Line : " << __LINE__ << ::std::endl;

#endif  //  HEU_MACROS_HPP