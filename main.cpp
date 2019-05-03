#include <QCoreApplication>
#include <QDebug>
#include "normaldist.h"
#include "rationalcubic.h"
#include "lbr.h"
#include<QDate>

//#include <QElapsedTimer>
//#define CONCAT_(x,y) x##y
//#define CONCAT(x,y) CONCAT_(x,y)

//#define CHECKTIME(x)  \
//    QElapsedTimer CONCAT(sb_, __LINE__); \
//    CONCAT(sb_, __LINE__).start(); \
//    x \
//    qDebug() << __FUNCTION__ << ":" << __LINE__ << " Elapsed time: " <<  CONCAT(sb_, __LINE__).elapsed() << " ms.";

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    int callput = -1; //call = 1; put = -1
    double strike = 200.0;
    double maturity = (20.0/365.0);
    double ttm = static_cast<double>(QDate::currentDate().daysTo(QDate(2019,5,24)))/365.0;
    qDebug()<<"calendar days to expiration:"<<static_cast<double>(QDate::currentDate().daysTo(QDate(2019,5,24)));
    double S = 204.3;
    double r = .025;
    double divPayabe = 0.73; // quarterly dividend
    double divTime = (12.0 / 365.0);
    double divYield = 0.0;
    double oPx = 4.45;
    double ivGuess = .30;
    double adjS = S * std::exp(-r * maturity);
    double divAdjS = adjS - divPayabe * std::exp(-r * divTime);

    //oPx,S,strike,maturity,callput[1,-1]
    lbr lbrObj;
    double iv = lbrObj.implied_volatility_from_a_transformed_rational_guess(oPx,S,strike,ttm,callput);
    qDebug() <<"r:"<<r<< "adjS:" << divAdjS << "iv:" << iv;
    return a.exec();
}
