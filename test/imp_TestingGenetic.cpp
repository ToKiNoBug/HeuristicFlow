/*
 Copyright © 2021  TokiNoBug
This file is part of OptimTemplates.

    OptimTemplates is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OptimTemplates is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OptimTemplates.  If not, see <https://www.gnu.org/licenses/>.

*/

#include "def_TestingGenetic.h"
#include <iostream>
#include <Eigen/Dense>
#include <ctime>
using namespace OptimT;
using namespace std;
/*
double randD(const double min,const double max) {
    return (max-min)*randD()+min;
}
*/
void testAckley_withRecord() {
    GAOption opt;
    opt.populationSize=200;
    opt.maxFailTimes=-1;
    opt.maxGenerations=3000;
    static const uint8_t MinIdx=0,MaxIdx=1,LrIdx=2;
    SOGA<array<double,2>,
            OptimT::FITNESS_LESS_BETTER,
            OptimT::RECORD_FITNESS,
            array<double,2>,//min
            array<double,2>,//max
            double//learning rate
            > algo;

    tuple<array<double,2>,array<double,2>,double> args;
    get<MinIdx>(args)={-5,-5};
    get<MaxIdx>(args)={5,5};
    get<LrIdx>(args)=0.05;

    algo.initialize(
                //initializer for array<double,2>
    [](array<double,2>* x,const typeof(args) * a) {
        for(uint32_t idx=0;idx<x->size();idx++) {
            x->operator[](idx)=randD(get<MinIdx>(*a)[idx],get<MaxIdx>(*a)[idx]);
        }},
    //Ackely function
    [](const array<double,2>* _x,const typeof(args) *,double * f) {
        double x=_x->operator[](0),y=_x->operator[](1);
        *f= -20*exp(-0.2*sqrt(0.5*(x*x+y*y)))
                -exp(0.5*(cos(M_2_PI*x)+cos(M_2_PI*y)))
                +20+M_E;},
    //crossover
    [](const array<double,2>* x,const array<double,2>* y,
            array<double,2> *X,array<double,2>*Y,
            const typeof(args) *) {
        const array<double,2> &copyx=*x,&copyy=*y;
        for(uint32_t idx=0;idx<x->size();idx++) {
            if(rand()%2)
                X->operator[](idx)=copyy[idx];
            else {
                X->operator[](idx)=copyx[idx];
            }
            if(rand()%2)
                Y->operator[](idx)=copyx[idx];
            else {
                Y->operator[](idx)=copyy[idx];
            }
        }},
    //mutate
    [](array<double,2>* x,const typeof(args) * a) {
        uint8_t mutateIdx=rand()%x->size();
        x->operator[](mutateIdx)+=get<LrIdx>(*a)*randD(-1,1);
        if(rand()%2)
            x->operator[](mutateIdx)*=-1;
        for(uint32_t idx=0;idx<x->size();idx++) {
            auto & var=x->operator[](idx);
            var=min(var,get<MaxIdx>(*a)[idx]);
            var=max(var,get<MinIdx>(*a)[idx]);
        }
    },
    //no other options
    nullptr,
    opt,
    args);

    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;

    cout<<"Solving spend "<<algo.generation()<<" generations in "
       <<double(t)/CLOCKS_PER_SEC<<" sec\n";
    cout<<"Result = ["<<algo.result()[0]<<" , "<<algo.result()[1]<<"]\n";

    cout<<"Fitness history :\n";

    for(auto i : algo.record()) {
        cout<<i<<'\n';
    }
    cout<<endl;

    //cout<<"Result idx = "<<algo.eliteIdx()<<endl;
}


void testSingleNumber() {
    SOGA<double,
            OptimT::FITNESS_GREATER_BETTER,
            OptimT::DONT_RECORD_FITNESS> algo;
    GAOption opt;
    opt.crossoverProb=0.8;
    opt.mutateProb=0.05;
    opt.maxFailTimes=100;
    opt.populationSize=100;
    opt.maxGenerations=3000;
    algo.initialize(
        [](double * x,const tuple<>*){*x=randD();},
        [](const double * x,const tuple<>*,double * f){*f=1.0/(abs(*x-3.0)+1e-8);},
        [](const double * x,const double *y,double *X,double *Y,const tuple<>*)
            {double mean=(*x+*y)/2;*X=(*x+mean)/2;*Y=(*y+mean)/2;},
        [](double * x,const tuple<>*)
            {
        //if(rand()%2)*x*=-1;
             *x+=randD(-1,1)*0.5;},
    nullptr,
    opt,
    tuple<>()
    );

    cout << "Start" << endl;
    algo.run();
    cout<<"Finished with "<<algo.generation()<<" generations"<<endl;
    cout<<"result = "<<algo.result()<<endl;
}

void testWithEigenLib() {
    Eigen::Array4d target(1,2,3,4);
    target/=target.sum();
    cout<<"Target="<<target.transpose()<<endl;
    GAOption opt;
    opt.populationSize=100;
    opt.maxGenerations=3000;
    opt.maxFailTimes=50;
    opt.crossoverProb=0.8;
    opt.mutateProb=0.05;

    SOGA<Eigen::Array4d,
            FITNESS_GREATER_BETTER,
            RECORD_FITNESS,
            Eigen::Array4d,Eigen::Array4d,Eigen::Array4d,double> algo;
    //val min max learning_rate
    tuple<Eigen::Array4d,Eigen::Array4d,Eigen::Array4d,double> Arg;
    static const uint8_t TargetOffset=0,MinOffset=1,MaxOffset=2,LROffset=3;
    get<TargetOffset>(Arg)=target;
    get<MinOffset>(Arg).setConstant(-1);
    get<MaxOffset>(Arg).setConstant(2);
    get<LROffset>(Arg)=0.01;


    algo.initialize(
                    [](Eigen::Array4d* x,const typeof(Arg)*){x->setRandom();},
    [](const Eigen::Array4d* x,const typeof(Arg)* arg,double *f){
        *f = -(*x-get<TargetOffset>(*arg)).square().maxCoeff();
    },
    [](const Eigen::Array4d*x,const Eigen::Array4d*y,
            Eigen::Array4d*X,Eigen::Array4d*Y,
            const typeof(Arg)*) {
        for(uint32_t i=0;i<4;i++) {
            X->operator()(i)=
                    (rand()%2)?
                        x->operator()(i):y->operator()(i);
            Y->operator()(i)=
                    (rand()%2)?
                        x->operator()(i):y->operator()(i);
        }},
    [](Eigen::Array4d*x,const typeof(Arg)* arg) {
        uint32_t idx=rand()%4;
        x->operator()(idx)+=randD(-1,1)*get<LROffset>(*arg);

        if(x->operator()(idx)>get<MaxOffset>(*arg)(idx)) {
            x->operator()(idx)=get<MaxOffset>(*arg)(idx);
        }
        if(x->operator()(idx)<get<MinOffset>(*arg)(idx)) {
            x->operator()(idx)=get<MinOffset>(*arg)(idx);
        }},
    nullptr,
    opt,
    Arg);


    algo.run();
    cout<<"Solving spend "<<algo.generation()<<" generations\n";
    cout<<"Result = "<<algo.result().transpose()<<endl;
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
            const vector<Point_t>*> algo;
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

    auto crossoverFun=[](const vector<double>*p1,const vector<double>*p2,
                      vector<double>*c1,vector<double>*c2,const Args_t*) {
        const uint32_t permL=p1->size();
        c1->resize(permL);c2->resize(permL);
        for(uint32_t i=0;i<permL;i++) {
            c1->operator[](i)=(std::rand()%2)?p1->operator[](i):p2->operator[](i);
            c2->operator[](i)=(std::rand()%2)?p1->operator[](i):p2->operator[](i);
        }
    };

    auto mutateFun=[](vector<double>*x,const Args_t*) {
        const uint32_t permL=x->size();
        if(randD()<0.5) {//modify one element's value
            double & v=x->operator[](std::rand()%permL);
            v+=randD(-0.01,0.01);
            v=min(v,1.0);
            v=max(v,0.0);
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
    algo.initialize(initializeFun,calculateFun,crossoverFun,mutateFun,nullptr,
                    opt,args);

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

    OptimT::NSGA2<std::array<double,XNum>,
            2,
            Std,
            FITNESS_LESS_BETTER,
            RECORD_FITNESS,
            PARETO_FRONT_DONT_MUTATE> algo;

    void (*iFun)(std::array<double,XNum>*,const std::tuple<>*) =
    [] (std::array<double,XNum> * x,const std::tuple<>*) {
        for(size_t i=0;i<XNum;i++) {
            x->operator[](i)=OptimT::randD();
        }
    };

    void (*fFun)
            (const std::array<double,XNum>* x,const std::tuple<>*, std::array<double,2> *f)
            =[](const std::array<double,XNum>* x,const std::tuple<>*, std::array<double,2> *f) {
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

    void (*cFun)(const std::array<double,XNum>*,const std::array<double,XNum>*,
                 std::array<double,XNum>*,std::array<double,XNum>*,const std::tuple<>*)
            =[](const std::array<double,XNum>*p1,const std::array<double,XNum>*p2,
            std::array<double,XNum>*ch1,std::array<double,XNum>*ch2,const std::tuple<>*)
    {
        for(size_t i=0;i<XNum;i++) {
            //discrete
            /*
            ch1->operator[](i)=(OptimT::randD()<0.5)?p1->operator[](i):p2->operator[](i);
            ch2->operator[](i)=(OptimT::randD()<0.5)?p1->operator[](i):p2->operator[](i);
            */
            ch1->operator[](i)=r*p1->operator[](i)+(1-r)*p2->operator[](i);
            ch2->operator[](i)=r*p2->operator[](i)+(1-r)*p1->operator[](i);
        }

    };

    void (*mFun)(std::array<double,XNum>*,const std::tuple<>*)=
            [](std::array<double,XNum>*x,const std::tuple<>*){
        const size_t mutateIdx=size_t(OptimT::randD(0,XNum))%XNum;

        x->operator[](mutateIdx)+=0.005*OptimT::randD(-1,1);

        x->operator[](mutateIdx)=std::min(x->operator[](mutateIdx),1.0);
        x->operator[](mutateIdx)=std::max(x->operator[](mutateIdx),0.0);
    };

    GAOption opt;
    opt.populationSize=2000;
    opt.maxFailTimes=500;
    opt.maxGenerations=2000;

    algo.initialize(iFun,fFun,cFun,mFun,nullptr,opt);
    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<"seconds and "<<algo.generation()<<"generations"<<endl;
    std::vector<std::array<double,2>> paretoFront;
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
            Std,
            FITNESS_LESS_BETTER,
            DONT_RECORD_FITNESS,
            PARETO_FRONT_DONT_MUTATE> algo;
    auto iFun=[](std::array<double,3> * x,const std::tuple<>*) {
        for(auto & i : *x) {
            i=randD(-5,5);
        }
    };
    auto fFun=[](const std::array<double,3> * x,const std::tuple<>*,std::array<double,2> *f) {
        double f1=0,f2=0;
        for(int i=0;i<2;i++) {
            f1+=-10*exp(-0.2*sqrt(OT_square(x->operator[](i))+OT_square(x->operator[](i+1))));
        }
        for(int i=0;i<3;i++) {
            f2+=pow(abs(x->operator[](i)),0.8)+5*sin(x->operator[](i)*x->operator[](i)*x->operator[](i));
        }
        f->operator[](0)=f1;
        f->operator[](1)=f2;
    };

    auto cFun=[](const std::array<double,3> *p1,const std::array<double,3> *p2,
            std::array<double,3> *ch1,std::array<double,3> *ch2,const std::tuple<>*) {
        for(int i=0;i<3;i++) {
            static const double r=0.2;
            ch1->operator[](i)=r*p1->operator[](i)+(1-r)*p2->operator[](i);
            ch2->operator[](i)=r*p2->operator[](i)+(1-r)*p1->operator[](i);
        }
    };

    auto mFun=[](std::array<double,3> * x,const std::tuple<>*) {
        const size_t idx=size_t(randD(0,3))%3;
        x->operator[](idx)+=0.1*randD(-1,1);
        x->operator[](idx)=min(x->operator[](idx),5.0);
        x->operator[](idx)=max(x->operator[](idx),-5.0);
    };

    GAOption opt;
    opt.maxGenerations=2000;
    opt.populationSize=1000;
    opt.maxFailTimes=-1;

    algo.initialize(iFun,fFun,cFun,mFun,nullptr,opt);

    algo.setCongestComposeFun(algo.default_ccFun_powered<3>);

    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<" seconds and "<<algo.generation()<<" generations"<<endl;
    std::vector<std::array<double,2>> paretoFront;
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
    using solver_t = 
    NSGA2<std::array<double,2>,
            2,
            DoubleVectorOption::Eigen,
            FITNESS_LESS_BETTER,
            DONT_RECORD_FITNESS,
            PARETO_FRONT_DONT_MUTATE>;
    solver_t algo;
    using Fitness_t = typename solver_t::Fitness_t;
    auto iFun=[](std::array<double,2> * x,const std::tuple<>*) {
        for(auto & i : *x) {
            i=randD(-5,5);
        }
    };
    auto fFun=[](const std::array<double,2> * _x,const std::tuple<>*,Fitness_t *f) {
        double & f1=f->operator[](0);
        double & f2=f->operator[](1);
        const double x=_x->operator[](0),y=_x->operator[](1);
        f1=4*(x*x+y*y);
        f2=OT_square(x-5)+OT_square(y-5);

        double constraint_g1=OT_square(x-5)+y*y-25;
        double constraint_g2=7.7-(OT_square(x-8)+OT_square(y+3));
        if(constraint_g1>0) {
            f1+=1e4+constraint_g1;
        }
        if(constraint_g2>0) {
            f2+=1e4+constraint_g2;
        }
    };

    auto cFun=[](const std::array<double,2> *p1,const std::array<double,2> *p2,
            std::array<double,2> *ch1,std::array<double,2> *ch2,const std::tuple<>*) {
        for(int i=0;i<2;i++) {
            static const double r=0.2;
            ch1->operator[](i)=r*p1->operator[](i)+(1-r)*p2->operator[](i);
            ch2->operator[](i)=r*p2->operator[](i)+(1-r)*p1->operator[](i);
        }
    };

    auto mFun=[](std::array<double,2> * x,const std::tuple<>*) {
        const size_t idx=size_t(randD(0,2))%2;
        x->operator[](idx)+=0.1*randD(-1,1);
        x->operator[](0)=min(x->operator[](0),5.0);
        x->operator[](0)=max(x->operator[](0),0.0);
        x->operator[](1)=min(x->operator[](1),3.0);
        x->operator[](1)=max(x->operator[](1),0.0);
    };

    GAOption opt;
    opt.maxGenerations=10000;
    opt.populationSize=200;
    opt.maxFailTimes=-1;

    algo.initialize(iFun,fFun,cFun,mFun,nullptr,opt);

    ///custom ccfun is not compulsory. Default value is nullptr
    //algo.setCongestComposeFun();

    cout<<"Start"<<endl;
    std::clock_t t=std::clock();
    algo.run();
    t=std::clock()-t;
    cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
       <<" seconds and "<<algo.generation()<<"generations"<<endl;

    /*
    std::vector<std::array<double,2>> paretoFront;
    algo.paretoFront(paretoFront);
    cout<<"paretoFront=[";
    for(const auto & i : paretoFront) {
        cout<<i[0]<<" , "<<i[1]<<";\n";
    }
    cout<<"];"<<endl;
    */

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
    Dynamic,
    DoubleVectorOption::Eigen,
    FITNESS_LESS_BETTER,
    DONT_RECORD_FITNESS,
    PARETO_FRONT_DONT_MUTATE,
    Eigen::Array<double,Dim,Eigen::Dynamic>>;

using Base_t = solver_t;
OptimT_MAKE_GABASE_TYPES
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

mutateFun mFun=[](Point_t * p,const ArgsType *) {
*p+=Point_t::Random()*0.05;
};

solver_t solver;
solver.setObjectiveNum(Dim);

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

solver.initialize(iFun,fFun,cFun,mFun,nullptr,opt,args);

clock_t t=clock();
solver.run();
t=std::clock()-t;
cout<<"Solving finished in "<<double(t)/CLOCKS_PER_SEC
    <<" seconds and "<<solver.generation()<<"generations"<<endl;

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