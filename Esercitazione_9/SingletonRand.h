#include "random.h"

class SingletonRand{
private:
    static SingletonRand * singleton_rnd;
    SingletonRand();
public:
    ~SingletonRand();
    Random * m_rnd;
    static SingletonRand * get_instance();
};

