/***
 *  $Id: iomanip.hpp 255 2009-06-15 14:38:15Z zola $
 **
 *  File: iomanip.hpp
 *  Created: Dec 12, 2007
 *
 *  Author: Jaroslaw Zola <jaroslaw.zola@gmail.com>
 */

#ifndef IOMANIP_HPP
#define IOMANIP_HPP

#include <mpix/MPI_env.hpp>
#include <iostream>


struct null_ostream : public std::ostream {
    null_ostream() : std::ios(0), std::ostream(0) { }
}; // struct null_ostream

struct mpi_iomanip {
    mpi_iomanip() : is_active(false), mpi_env(0) { }
    bool is_active;
    const mpix::MPI_env* mpi_env;
}; // struct AppIO

extern null_ostream null_out;
extern mpi_iomanip mpi_iom;


inline std::ostream& master(std::ostream& os) {
    return mpi_iom.mpi_env->am_I_root() ? os : null_out;
} // master

inline std::ostream& active(std::ostream& os) {
    return mpi_iom.is_active ? os : null_out;
} // active

inline std::ostream& every(std::ostream& os) {
    os << "[" << mpi_iom.mpi_env->rank() << "] ";
    return os;
} // every

class host {
public:
    explicit host(int id) : id_(id) { }

private:
    int id_;

    friend std::ostream& operator<<(std::ostream& os, const host& h) {
	if (mpi_iom.mpi_env->rank() != h.id_) return null_out;
	os << "[" << h.id_ << "] ";
	return os;
    } // operator <<

}; // class host

inline std::ostream& warning(std::ostream& os) {
    os << "Warning: ";
    return os;
} // warning

inline std::ostream& error(std::ostream& os) {
    os << "Error: ";
    return os;
} // error

inline std::ostream& oops(std::ostream& os) {
    os << "Oops!!! ";
    return os;
} // oops

class timer {
public:
    explicit timer(double t) : t_(t) { }

private:
    double t_;

    friend std::ostream& operator<<(std::ostream& os, const timer& t) {
	unsigned int tt = static_cast<unsigned int>(t.t_);
	unsigned int ht = tt / 3600;
	unsigned int mt = (tt % 3600) / 60;
	unsigned int st = (tt % 3600) % 60;
	os << t.t_ << "s (" << ht << "h" << mt << "m" << st << "s)";
	return os;
    } // operator <<

}; // class timer

#endif // IOMANIP_HPP
