
#ifndef GPUPOLY_EXPORT_H
#define GPUPOLY_EXPORT_H

#ifdef GPUPOLY_STATIC_DEFINE
#define GPUPOLY_EXPORT
#define GPUPOLY_NO_EXPORT
#else
#ifndef GPUPOLY_EXPORT
#ifdef gpupoly_EXPORTS
/* We are building this library */
#define GPUPOLY_EXPORT __attribute__((visibility("default")))
#else
/* We are using this library */
#define GPUPOLY_EXPORT __attribute__((visibility("default")))
#endif
#endif

#ifndef GPUPOLY_NO_EXPORT
#define GPUPOLY_NO_EXPORT __attribute__((visibility("hidden")))
#endif
#endif

#ifndef GPUPOLY_DEPRECATED
#define GPUPOLY_DEPRECATED __attribute__((__deprecated__))
#endif

#ifndef GPUPOLY_DEPRECATED_EXPORT
#define GPUPOLY_DEPRECATED_EXPORT GPUPOLY_EXPORT GPUPOLY_DEPRECATED
#endif

#ifndef GPUPOLY_DEPRECATED_NO_EXPORT
#define GPUPOLY_DEPRECATED_NO_EXPORT GPUPOLY_NO_EXPORT GPUPOLY_DEPRECATED
#endif

#if 0 /* DEFINE_NO_DEPRECATED */
#ifndef GPUPOLY_NO_DEPRECATED
#define GPUPOLY_NO_DEPRECATED
#endif
#endif

#endif /* GPUPOLY_EXPORT_H */
