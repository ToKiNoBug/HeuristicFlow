// This file is part of Eigen, a lightweight C++ template library
// for linear algebra.
//
// Copyright (C) 2022 Shawn Li <tokinobug@163.com>
//
// This Source Code Form is subject to the terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "def_TestingGenetic.h"
#include <iostream>
#include <Eigen/Dense>
#include <ctime>
using namespace Heu;
using namespace std;

void testAckley_withRecord() {

    using args_t = Heu::BoxNdS<2,Std>;

    using solver_t = 
    SOGA<array<double,2>,
            Heu::FITNESS_LESS_BETTER,
            Heu::RECORD_FITNESS,
            args_t,
    Heu::GADefaults<array<double,2>,Std,args_t>::iFunNd<>,
    nullptr,
    Heu::GADefaults<array<double,2>,Std,args_t>::cFunNd,
    Heu::GADefaults<array<double,2>,Std,args_t>::mFun_d<>>;
    solver_t algo;
    
    GAOption opt;
    opt.populationSize=50;
    opt.maxFailTimes=-1;
    opt.maxGenerations=100;
    
    algo.setOption(opt);

    {
        args_t args;
        args.setMin(-5);
        args.setMax(5);
        args.setLearnRate(0.05);
        algo.setArgs(args);
    }

    algo.setfFun(
    //Ackely function
    [](const array<double,2>* _x,const solver_t::ArgsType *,double * f) {
        double x=_x->operator[](0),y=_x->operator[](1);
        *f= -20*exp(-0.2*sqrt(0.5*(x*x+y*y)))
                -exp(0.5*(cos(M_2_PI*x)+cos(M_2_PI*y)))
                +20+M_E;}
    );

    algo.initializePop();

    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<algo.bestFitness();
/*
    cout<<"Solving spend "<<algo.generation()<<" generations in "
       <<double(t)/CLOCKS_PER_SEC<<" sec\n";
    cout<<"Result = ["<<algo.result()[0]<<" , "<<algo.result()[1]<<"]\n";
*/
    /*
    cout<<"Fitness history :\n";

    for(auto i : algo.record()) {
        cout<<i<<'\n';
    }
    cout<<endl;
    */
}


void testTSP(const uint32_t PointNum) {
    static const uint8_t DIM=2;
    //static const double LengthBase=100;
    typedef array<double,DIM> Point_t;
    vector<Point_t> points;
    points.clear();
    points.reserve(PointNum);
    for(uint32_t i=0;i<PointNum;i++) {
        points.emplace_back();
        for(auto & i : points.back()) {
            i=randD();
        }
    }
    /*
    cout<<"Generated random points :"<<endl;
    for(const auto & i : points) {
        cout<<'(';
        for(auto j : i) {
            cout<<j<<" , ";
        }
        cout<<")"<<endl;
    }
    */
    typedef pair<double,uint32_t> permUnit;

    //typedef vector<permUnit> permulation;
    //      var,        less=better,    data src
    SOGA<vector<double>,
            FITNESS_LESS_BETTER,
            DONT_RECORD_FITNESS,
            std::tuple<const vector<Point_t>*>> algo;
    static const uint8_t dataIdx=0;
    typedef tuple<const vector<Point_t>*> Args_t;
    Args_t args;//=make_tuple(PointNum,points.data());
    get<dataIdx>(args)=&points;
    //initialize function

    auto initializeFun=[](vector<double> * x,const Args_t* args) {
        const uint32_t permL=get<dataIdx>(*args)->size();
        x->resize(permL);
        for(uint32_t i=0;i<permL;i++) {
            x->operator[](i)=randD();
        }
    };

    //calculate fitness
    auto calculateFun=[](const vector<double> *x,const Args_t* args,double *f) {

        const uint32_t permL=x->size();
        vector<permUnit> perm(permL);
        for(uint32_t i=0;i<permL;i++) {
            perm[i].first=x->operator[](i);
            perm[i].second=i;
        }
        std::sort(perm.begin(),perm.end(),
                  //compare function for std::sort
                  [](const permUnit &a,const permUnit &b)
                      { return a.first>b.first;});

        double L=0;
        for(uint32_t i=1;i<permL;i++) {
            const Point_t &prev=get<dataIdx>(*args)->operator[](perm[i-1].second);
            const Point_t &cur=get<dataIdx>(*args)->operator[](perm[i].second);
            double curL=0;
            for(uint8_t d=0;d<DIM;d++) {
                curL+=(prev[d]-cur[d])*(prev[d]-cur[d]);
            }
            L+=curL;
        } *f= L;
    };

    //calculate natural perm pathL
    {
        vector<double> natural(PointNum);
        for(uint32_t i=0;i<PointNum;i++) {
            natural[i]=double(i)/PointNum;
        }
        double naturalPathL;
        calculateFun(&natural,&args,&naturalPathL);
        cout<<"default pathL = "<<naturalPathL<<endl;
    }

    auto crossoverFun
            =Heu::GADefaults<vector<double>,DoubleVectorOption::Std,Args_t>::
                cFunXd<Heu::DivCode::Half>;

    auto mutateFun=[](const vector<double>*src,vector<double>*x,const Args_t*) {
        *x=*src;
        const uint32_t permL=x->size();
        if(randD()<0.5) {//modify one element's value
            double & v=x->operator[](std::rand()%permL);
            v+=randD(-0.01,0.01);
            v=std::min(v,1.0);
            v=std::max(v,0.0);
        }
        else {
            const uint32_t flipB=std::rand()%(permL-1);
            const uint32_t flipE=std::rand()%(permL-flipB-1)+flipB+1;
            //cout<<"flipB="<<flipB<<" , flipE"<<flipE<<endl;
            const uint32_t flipN=flipE-flipB+1;
            for(uint32_t flipped=0;flipped*2<=flipN;flipped++) {
                swap(x->operator[](flipB+flipped),x->operator[](flipE-flipped));
            }
        }
    };


    GAOption opt;
    opt.maxGenerations=30*PointNum;
    opt.maxFailTimes=-1;

    algo.setiFun(initializeFun);
    algo.setmFun(mutateFun);
    algo.setfFun(calculateFun);
    algo.setcFun(crossoverFun);
    algo.setOption(opt);
    algo.setArgs(args);
    algo.initializePop();

    cout<<"run!"<<endl;
    std::clock_t c=std::clock();
    algo.run();
    c=std::clock()-c;
    cout<<"finished with "<<algo.generation()<<" generations and "
       <<double(c)/CLOCKS_PER_SEC<<" s\n";
    cout<<"result fitness = "<<algo.bestFitness()<<endl;
}

///Zitzler–Deb–Thiele's function N. 3
void testNSGA2_ZDT3() {
    //0<=x_i<=1, 1<=i<=30
    static const size_t XNum=30;
    static const double r=0.2;

    Heu::NSGA2<std::array<double,XNum>,
            2,
            FITNESS_LESS_BETTER,
            RECORD_FITNESS> algo;

    void (*iFun)(std::array<double,XNum>*) =
    [] (std::array<double,XNum> * x) {
        for(size_t i=0;i<XNum;i++) {
            x->operator[](i)=Heu::randD();
        }
    };

    void (*fFun)
            (const std::array<double,XNum>* x, Eigen::Array<double,2,1> *f)
            =[](const std::array<double,XNum>* x, Eigen::Array<double,2,1> *f) {
      f->operator[](0)=x->operator[](0);
      const double && f1=std::move(f->operator[](0));
      double g=0;
      for(size_t i=1;i<XNum;i++) {
          g+=x->operator[](i);
      }
      g=1+g*9.0/(XNum-1);
      double && f1_div_g=f1/g;
      double && h=1-std::sqrt(f1_div_g)-f1_div_g*std::sin(10*M_PI*f1);
      f->operator[](1)=g*h;
    };

    auto cFun=
            Heu::GADefaults<std::array<double,XNum>>::
                cFunNd<(Heu::encode<1,2>::code)>;

    void (*mFun)(const std::array<double,XNum>*src,std::array<double,XNum>*)=
            [](const std::array<double,XNum>*src,std::array<double,XNum>*x){
        *x=*src;
        const size_t mutateIdx=Heu::randIdx(XNum);

        x->operator[](mutateIdx)+=0.005*Heu::randD(-1,1);

        x->operator[](mutateIdx)=std::min(x->operator[](mutateIdx),1.0);
        x->operator[](mutateIdx)=std::max(x->operator[](mutateIdx),0.0);
    };

    GAOption opt;
    opt.populationSize=2000;
    opt.maxFailTimes=500;
    opt.maxGenerations=2000;

    algo.setiFun(iFun);
    algo.setmFun(mFun);
    algo.setfFun(fFun);
    algo.setcFun(cFun);
    algo.setOption(opt);
    algo.initializePop();

    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<"seconds and "<<algo.generation()<<"generations"<<endl;
    std::vector<Eigen::Array<double,2,1>> paretoFront;
    algo.paretoFront(paretoFront);
    cout<<"paretoFront=[";
    for(const auto & i : paretoFront) {
        cout<<i[0]<<" , "<<i[1]<<";\n";
    }
    cout<<"];"<<endl;
    /*
    cout<<"\n\n\n population=[";
    for(const auto & i : algo.population()) {
        cout<<i.fitness()[0]<<" , "<<i.fitness()[1]<<";\n";
    }
    cout<<"];"<<endl;
    */
}

void testNSGA2_Kursawe() {
    //1<=i<=3,  -5<=x_i<=5
    NSGA2<std::array<double,3>,
            2,
            FITNESS_LESS_BETTER,
            DONT_RECORD_FITNESS> algo;
    auto iFun=[](std::array<double,3> * x) {
        for(auto & i : *x) {
            i=randD(-5,5);
        }
    };
    auto fFun=[](const std::array<double,3> * x,Eigen::Array<double,2,1> *f) {
        double f1=0,f2=0;
        for(int i=0;i<2;i++) {
            f1+=-10*exp(-0.2*sqrt(square(x->operator[](i))+square(x->operator[](i+1))));
        }
        for(int i=0;i<3;i++) {
            f2+=pow(abs(x->operator[](i)),0.8)+5*sin(x->operator[](i)*x->operator[](i)*x->operator[](i));
        }
        f->operator[](0)=f1;
        f->operator[](1)=f2;
    };

    auto cFun=[](const std::array<double,3> *p1,const std::array<double,3> *p2,
            std::array<double,3> *ch1,std::array<double,3> *ch2) {
        for(int i=0;i<3;i++) {
            static const double r=0.2;
            ch1->operator[](i)=r*p1->operator[](i)+(1-r)*p2->operator[](i);
            ch2->operator[](i)=r*p2->operator[](i)+(1-r)*p1->operator[](i);
        }
    };

    auto mFun=[](const std::array<double,3> * src,std::array<double,3> * x) {
        *x=*src;
        const size_t idx=randIdx(3);
        x->operator[](idx)+=0.1*randD(-1,1);
        x->operator[](idx)=std::min(x->operator[](idx),5.0);
        x->operator[](idx)=std::max(x->operator[](idx),-5.0);
    };

    GAOption opt;
    opt.maxGenerations=2000;
    opt.populationSize=1000;
    opt.maxFailTimes=-1;

    algo.setiFun(iFun);
    algo.setmFun(mFun);
    algo.setfFun(fFun);
    algo.setcFun(cFun);
    algo.setOption(opt);
    algo.initializePop();

    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<" seconds and "<<algo.generation()<<" generations"<<endl;
    std::vector<Eigen::Array<double,2,1>> paretoFront;
    algo.paretoFront(paretoFront);
    cout<<"paretoFront=[";
    for(const auto & i : paretoFront) {
        cout<<i[0]<<" , "<<i[1]<<";\n";
    }
    cout<<"];"<<endl;
    /*
    cout<<"\n\n\n population=[";
    for(const auto & i : algo.population()) {
        cout<<i.fitness()[0]<<" , "<<i.fitness()[1]<<";\n";
    }
    cout<<"];"<<endl;
    */
}

void testNSGA2_Binh_and_Korn() {
    //0<=x_0<=5,  0<=x_1<=3

    using args_t = Heu::BoxNdN<2,Heu::DoubleVectorOption::Std>;

    using solver_t = 
    NSGA2<std::array<double,2>,
            2,
            FITNESS_LESS_BETTER,
            RecordOption::DONT_RECORD_FITNESS,args_t,
            Heu::GADefaults<std::array<double,2>,Std,args_t>::iFunNd,
            nullptr,
            Heu::GADefaults<std::array<double,2>,Std,args_t>::cFunNd<>,
            Heu::GADefaults<std::array<double,2>,Std,args_t>::mFun_d            
            >;

    solver_t algo;

    using Fitness_t = typename solver_t::Fitness_t;

    auto fFun=[](const std::array<double,2> * _x,const args_t *,Fitness_t *f) {
        double f1;
        double f2;
        const double x=_x->operator[](0),y=_x->operator[](1);
        f1=4*(x*x+y*y);
        f2=square(x-5)+square(y-5);

        double constraint_g1=square(x-5)+y*y-25;
        double constraint_g2=7.7-(square(x-8)+square(y+3));

        if(constraint_g1>0) {
            f1=1e4+constraint_g1;
        }
        if(constraint_g2>0) {
            f2=1e4+constraint_g2;
        }

        *f={f1,f2};
    };

    {
        GAOption opt;
        opt.maxGenerations=400;
        opt.populationSize=200;
        opt.maxFailTimes=-1;
        opt.crossoverProb=0.8;
        opt.mutateProb=0.1;
        algo.setOption(opt);
    }
    {
        args_t box;
        box.setMin({0,0});
        box.setMax({5,3});
        box.setLearnRate({0.05,0.03});
        algo.setArgs(box);
    }

    algo.setfFun(fFun);
    algo.initializePop();

    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<" seconds and "<<algo.generation()<<"generations"<<endl;

    
    cout<<"paretoFront=[";
    for(const auto & i : algo.pfGenes()) {
        cout<<i->_Fitness[0]<<" , "<<i->_Fitness[1]<<";\n";
    }
    cout<<"];"<<endl;
    

    /*
    cout<<"\n\n\n population=[";
    for(const auto & i : algo.population()) {
        cout<<i.fitness()[0]<<" , "<<i.fitness()[1]<<";\n";
    }
    cout<<"];"<<endl;
    */
}


void testDistance(const size_t N) {
static const size_t Dim=2;
using Point_t = Eigen::Array<double,Dim,1>;
using solver_t = NSGA2<Point_t,
    Runtime,
    FITNESS_LESS_BETTER,
    DONT_RECORD_FITNESS,
    std::tuple<Eigen::Array<double,Dim,Eigen::Dynamic>>>;

using Base_t = solver_t;
Heu_MAKE_NSGABASE_TYPES
using Fitness_t = solver_t::Fitness_t;

static const double r=0.2;

initializeFun iFun=[](Point_t *p,const ArgsType *) {
p->setRandom();
};

fitnessFun fFun=[](const Point_t * p,const ArgsType* arg,Fitness_t *f) {
const size_t dotNum=std::get<0>(*arg).cols();
auto diff=(std::get<0>(*arg)-(p->replicate(1,dotNum))).square();

*f=diff.colwise().sum();

};

crossoverFun cFun=[](const Point_t * p1,const Point_t * p2,
    Point_t * c1,Point_t * c2,
    const ArgsType*) {
    *c1=r*(*p1)+(1-r)*(*p2);
    *c2=r*(*p2)+(1-r)*(*p1);
};

mutateFun mFun=[](const Point_t *src,Point_t * p,const ArgsType *) {
    *p=*src;
    *p+=Point_t::Random()*0.05;
};

solver_t solver;
solver.setObjectiveNum(N);

GAOption opt;
opt.maxGenerations=10000;
opt.maxFailTimes=4000;
opt.populationSize=400;

ArgsType args;
std::get<0>(args)=(Eigen::Array<double,Dim,Eigen::Dynamic>::Random(Dim,N))/2+0.5;

for(size_t c=1;c<N;c++) {
    auto prevC=std::get<0>(args).col(c-1);
    std::get<0>(args).col(c)=4*prevC*(1-prevC);
}

std::get<0>(args)=(std::get<0>(args)-0.5);

cout<<"args=\n";
cout<<std::get<0>(args)<<endl;

solver.setiFun(iFun);
solver.setcFun(cFun);
solver.setfFun(fFun);
solver.setmFun(mFun);
solver.setArgs(args);
solver.setOption(opt);

solver.initializePop();

clock_t t=clock();
solver.run();
t=std::clock()-t;
cout<<"\n\nSolving finished in "<<double(t)/CLOCKS_PER_SEC
    <<" seconds and "<<solver.generation()<<" generations\n\n"<<endl;

cout<<"PFValue=[";
for(const auto & i : solver.pfGenes()) {
    cout<<i->_Fitness.transpose()<<'\n';
}
cout<<"];"<<endl;

cout<<"\n\n\n\nPFPos=[";
for(const auto & i : solver.pfGenes()) {
    cout<<i->self.transpose()<<'\n';
}
cout<<"];"<<endl;

}

template<size_t M,size_t N,size_t realM=M>
void DTLZ7(const Eigen::Array<double,N,1> * x,EigenVecD_t<M> * f) {
    static_assert(realM>=2,"actural objective amount mustn't be less than 2");
    f->resize(realM,1);
    auto xm=x->template bottomRows<N-realM+1>();
    f->template topRows<realM-1>()=x->template topRows<realM-1>();
    const double g=1+9*(xm.sum())/
            //std::sqrt(1e-40+xm.square().sum())
            (N-realM+1)
            ;
    auto fi=f->template topRows<realM-1>();
    auto one_add_sum3pifi=1+(3*M_PI*fi).sin();
    const double h=realM-(fi*one_add_sum3pifi/(1+g)).sum();
    f->operator[](realM-1)=(1+g)*h;
};


void testNSGA3_DTLZ7() {
static const size_t N=20;
static const size_t M=6;
using solver_t = NSGA3<Eigen::Array<double,N,1>,
    M,
    //FITNESS_LESS_BETTER,
    DONT_RECORD_FITNESS,
    DOUBLE_LAYER,
    void>;

using Var_t = Eigen::Array<double,N,1>;
using Fitness_t = solver_t::Fitness_t;

auto iFun=GADefaults<Var_t,DoubleVectorOption::Eigen,void>::iFunNd<>;

auto DTLZ1=[](const Var_t * x,Fitness_t * f) {
    f->resize(M,1);
    auto xm=x->bottomRows<N-M+1>();
    auto xm_sub_half_square=(xm-0.5).square();
    auto cos_20pi_mul=(20*M_PI*(xm-0.5)).cos();
    const double g=100*(std::sqrt(xm.square().sum())
        +(xm_sub_half_square-cos_20pi_mul).sum());
    double accum=0.5*(1+g);
    int64_t dim=0;
    for(int64_t obj=M-1;obj>=0;obj--) {
        if(obj>0) {
            (*f)[obj]=accum*(1-(*x)[dim]);
            accum*=(*x)[dim];
            dim++;
        }
        else {
            (*f)[obj]=accum;
        }
    }
};


auto cFun=GADefaults<Var_t,DoubleVectorOption::Eigen>::cFunNd<encode<1,10>::code>;


auto mFun=[](const Var_t * src,Var_t * v) {
    *v=*src;
    const size_t idx=randIdx(v->size());
    double & p=v->operator[](idx);
    p+=0.5*randD(-1,1);
    p=std::min(p,1.0);
    p=std::max(p,0.0);
};

GAOption opt;
cout<<"maxGenerations=";
cin>>opt.maxGenerations;

opt.maxFailTimes=-1;
cout<<"populationSize=";
cin>>opt.populationSize;
opt.crossoverProb=0.8;
opt.mutateProb=0.1;

solver_t solver;
//solver.setObjectiveNum(M);
solver.setiFun(iFun);
solver.setfFun(DTLZ7<M,N,M>);
solver.setcFun(cFun);
solver.setmFun(mFun);
solver.setOption(opt);
solver.setReferencePointPrecision(6,4); cout<<"RPCount="<<solver.referencePointCount()<<endl;
solver.initializePop();

//cout<<"RP=["<<solver.referencePoints()<<"]';\n\n\n"<<endl;


clock_t c=clock();
solver.run();
c=clock()-c;

cout<<"solving finished in "<<c<<" ms with "<<solver.generation()<<" generations."<<endl;

cout<<"PFV=[";
for(const auto & i : solver.pfGenes()) {
    cout<<i->_Fitness.transpose()<<";\n";
}
cout<<"];\n\n\n"<<endl;

}

template<size_t M,size_t N,size_t precision>
void searchPFfun(
    size_t nIdx,
    Eigen::Array<double,N,1> * var,
    vector<pair<Eigen::Array<double,M,1>,size_t>> * dst) {

    static_assert(precision>=2,"You should assign at least 2 on a single dim");
    if(nIdx+1>=N) {
        Eigen::Array<double,M,1> f;
        for(size_t i=0;i<=precision;i++) {
            var->operator[](nIdx)=double(i)/precision;
            DTLZ7<M,N>(var,&f);
            dst->emplace_back(make_pair(f,0ULL));
        }
        return;
    }

    for(size_t i=0;i<=precision;i++) {
        var->operator[](nIdx)=double(i)/precision;
        searchPFfun<M,N,precision>(nIdx+1,var,dst);
    }

}

void searchPF() {
    static const size_t N=5,M=3;
    auto fun=DTLZ7<M,N>;
    using Var_t = Eigen::Array<double,N,1>;
    using Fitness_t = Eigen::Array<double,M,1>;
    using pair_t = pair<Fitness_t,size_t>;
    vector<pair_t> pop;
    static const size_t precision=10;
    pop.reserve(size_t(pow(N,precision)));
    Var_t v;

    cout<<"traversing"<<endl;
    clock_t c=clock();
    searchPFfun<M,N,precision>(0,&v,&pop);
    c=clock()-c;
    cout<<"traversing cost "<<c<<" ms"<<endl;

    cout<<"Nondominated sorting"<<endl;
    c=clock();
    static const int32_t thN=HfGlobal::threadNum();
#pragma omp parallel for
    for(int32_t begIdx=0;begIdx<thN;begIdx++) {
        for(int32_t idx=begIdx;idx<pop.size();idx+=thN) {
            pop[idx].second=0;
            for(size_t er=0;er<pop.size();er++) {
                pop[idx].second
                    +=internal::Pareto<M,FITNESS_LESS_BETTER>::
                    isStrongDominate(&pop[er].first,&pop[idx].first);
            }
        }
    }
    c=clock()-c;
    cout<<"NS cost "<<c<<" ms"<<endl;

    cout<<"StdPFV=[";
    for(const auto & i : pop) {
        if(i.second==0)
            cout<<i.first.transpose()<<";\n";
    }
    cout<<"];\n\n\n"<<endl;
}
