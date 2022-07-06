#ifndef HEU_TESTFUNCTIONSCOMMON_HPP
#define HEU_TESTFUNCTIONSCOMMON_HPP
namespace heu {
namespace internal {

#define HEU_REPEAT_FUNCTIONS(className, functionName)                                     \
  inline static void functionName(const Var_t *x, const Arg_t *, Fitness_t *f) noexcept { \
    className<Var_t, Fitness_t, void>::functionName(x, f);                                \
  }

}  // namespace internal
}  // namespace heu
#endif  // HEU_TESTFUNCTIONSCOMMON_HPP