#include "mainwindow.h"
#include "ui_mainwindow.h"

static std::random_device rd;           //Генератор истинно случайных чисел (шум из аппаратных устр-в)
static std::default_random_engine gen;  //Генератор псевдорандомных чисел


static QVector<double> val;

void MainWindow::paintGraph(QCustomPlot* dst,QVector<double>& X,QVector<double>& F1,QString name1,
                            QVector<double>& Xaux, QVector<double>& F2, QString name2)
{
    double minF = F1[0], maxF = F1[0];
    for(int i = 1; i < F1.size(); i++)
    {
        if(F1[i]<minF) minF = F1[i];
        if(F1[i]>maxF) maxF = F1[i];
    }

    dst->clearPlottables();

    QPen pen;
    pen.setColor(Qt::GlobalColor(7));

    QCPCurve* F_sample;
    F_sample = new QCPCurve(dst->xAxis, dst->yAxis);
    F_sample->setData(Xaux,F2);
    F_sample->setPen(pen);

    dst->addGraph();
    dst->graph(0)->setData(X,F1);
    dst->graph(0)->setName(QString(name1));
    dst->plottable(0)->setName(QString(name2));
    dst->legend->setVisible(true);
    dst->xAxis->setRange(X[0],X[X.size()-1]);
    dst->yAxis->setRange(minF,maxF);
    dst->replot();
}

//Интегральная функция распределения
double intF(double lambda, double x){
    if(x>=0)
        return 1-exp(-lambda*x);
    else
        return 0;
}

//Функция плотности
double densf(double lambda, double x){
    if(x>=0)
        return lambda*exp(-lambda*x);
    else
        return 0;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow) {
    ui->setupUi(this);
    on_histN_textChanged(ui->histN->text());
    on_hypN_textChanged(ui->hypN->text());
    ui->rvTable->setColumnWidth(0,50);
    ui->rvTable->horizontalHeader()->setVisible(false);
    ui->rvStatsTable->verticalHeader()->setVisible(false);
    ui->histRanges->horizontalHeader()->setVisible(false);
    ui->histStats->horizontalHeader()->setVisible(false);
    ui->hypRanges->horizontalHeader()->setVisible(false);
    ui->hypH0->horizontalHeader()->setVisible(false);
    ui->distrFunc->setInteraction(QCP::iRangeZoom,true);
    ui->distrFunc->setInteraction(QCP::iRangeDrag, true);
    ui->histGraph->setInteraction(QCP::iRangeZoom,true);
    ui->histGraph->setInteraction(QCP::iRangeDrag, true);
}

MainWindow::~MainWindow() {
    delete ui;
}

void MainWindow::on_startButton_clicked() {
    ui->Params->setDisabled(true);
    QTableWidgetItem* tbl;
    QVector<double> unX,valX;
    QVector<double> F,F_sample;
    unsigned int seed;

    //подготовка для рандомайзера
    int maxrand = ui->uniqValN->text().toInt();
    if(ui->isrdSeed->checkState() == 2){
        seed = rd();
        ui->curSeed->setText(QString::number(seed));
    } else
        seed = ui->curSeed->text().toUInt();
    gen.seed(seed);
    std::uniform_int_distribution<> dist(1,maxrand);

    double M,Avg = 0,M_Avg,D,S2 = 0,D_S2,Me,R;
    double lambda = ui->lambdaPar->text().toDouble();
    int n = ui->nPar->text().toInt();


    val.resize(n);
    valX.resize(2*n);
    F_sample.resize(2*n);
    ui->rvTable->setColumnCount(n);

    //Получение n значений с.в. с экспоненц. распред.
    for(int i = 0;i<n;i++) {
        val[i] = -log((double)dist(gen)/maxrand)/lambda;
        Avg+=val[i];
    }

    std::sort(val.begin(),val.end());   //Упорядочивание значений с.в. по возрастанию
    M = 1/lambda;                       //Матожидание экспоненц. распред.
    Avg/=n;                             //Выборочное среднее
    M_Avg = std::abs(M-Avg);                 //Модуль разности M и Выб. сред.
    D = 1/(lambda*lambda);              //Дисперсия экспоненц. распред.
    for(int i = 0;i<n;i++)
        S2+=(val[i]-Avg)*(val[i]-Avg);
    S2/=n;                              //Выборочная дисперсия
    D_S2 = std::abs(D-S2);                   //Модуль разности D и Выб. дисп.
    if(n%2 == 0)                        //Выб. медиана
        Me = (val[n/2-1]+val[n/2])/2;
    else
        Me = val[n/2];
    R = val[n-1]-val[0];                //Размах выборки

    //Вычисление функц. распр. и выборочной функц. распр.

    //для построения функции распределения сетка n+1000 чтобы красиво
    for(double arg = val[0];arg<=val[n-1];arg+=((val[n-1]-val[0])/(n+1000))){
        unX.push_back(arg);
        F.push_back(intF(lambda,arg));
    }

    //вычислить выборочную функцию и добавить фейковые точки для ступенчатости
    for(int i = 0; i<n;i++) {
        valX[2*i] = val[i];
        if(i!=n-1)
            valX[2*i+1] = val[i+1];
        else
            valX[2*i+1] = val[n-1]+0.1;
        F_sample[2*i] = F_sample[2*i+1] = (double)(i+1)/n;
    }

    //Вычисление меры расхождения
    double maxDev = 0;
    for(int i = 0; i<n;i++) {
        double tempmaxDev;
        if(i!=n-1)
            tempmaxDev = std::max(std::abs(intF(lambda,val[i])-F_sample[2*i]),std::abs(intF(lambda,val[i+1]) - F_sample[2*i+1]));
        else
            tempmaxDev = std::max(std::abs(F[i]-F_sample[2*i]),0.0);
        if(tempmaxDev>maxDev)
            maxDev = tempmaxDev;
    }

    //отрисовать графики функций распределения
    paintGraph(ui->distrFunc,unX,F,"F(x)",valX,F_sample,"F выбор.");

    //заполнить таблицу характеристик с.в.
    tbl = new QTableWidgetItem(QString::number(M));
    ui->rvStatsTable->setItem(0,0,tbl);
    tbl = new QTableWidgetItem(QString::number(Avg));
    ui->rvStatsTable->setItem(0,1,tbl);
    tbl = new QTableWidgetItem(QString::number(M_Avg));
    ui->rvStatsTable->setItem(0,2,tbl);
    tbl = new QTableWidgetItem(QString::number(D));
    ui->rvStatsTable->setItem(0,3,tbl);
    tbl = new QTableWidgetItem(QString::number(S2));
    ui->rvStatsTable->setItem(0,4,tbl);
    tbl = new QTableWidgetItem(QString::number(D_S2));
    ui->rvStatsTable->setItem(0,5,tbl);
    tbl = new QTableWidgetItem(QString::number(Me));
    ui->rvStatsTable->setItem(0,6,tbl);
    tbl = new QTableWidgetItem(QString::number(R));
    ui->rvStatsTable->setItem(0,7,tbl);
    tbl = new QTableWidgetItem(QString::number(maxDev));
    ui->rvStatsTable->setItem(0,8,tbl);
    for(int i = 0;i<n;i++){
        tbl = new QTableWidgetItem(QString::number(val[i]));
        ui->rvTable->setItem(0,i,tbl);
    }

    //поправить раздел гистограммы
    on_histN_textChanged(ui->histN->text());
    on_hypN_textChanged(ui->hypN->text());
    ui->max->setText("");
    on_histStart_clicked();
    on_hypStart_clicked();
    ui->Params->setDisabled(false);
}

//При изменении выбора сида
void MainWindow::on_isrdSeed_stateChanged(int state)
{
    if(state == 0)
        ui->curSeed->setDisabled(false);
    else
         ui->curSeed->setDisabled(true);
}



//Построение гистограммы
void MainWindow::on_histStart_clicked()
{
    if(!val.isEmpty()){

        double lambda = ui->lambdaPar->text().toDouble();
        double maxdif = 0;
        int n =val.size();
        QTableWidgetItem* tbl;
        QVector<double> unX,valX;
        QVector<double> f, hist;

    //для построения функции плотности сетка n+1000 чтобы красиво
        for(double arg = val[0];arg<=val[n-1];arg+=((val[n-1]-val[0])/(n+1000))){
            unX.push_back(arg);
            f.push_back(densf(lambda,arg));
        }

        QVector<double> ranges,nums;
        //Вытащить границы разбиения из таблицы
        ranges.push_back(val[0]);
        for(int i = 1; i < ui->histRanges->columnCount()-1; i++)
            ranges.push_back(ui->histRanges->item(0,i)->text().toDouble());
        ranges.push_back(val[val.size()-1]);
        int k = 0;
        nums.resize(ranges.size()-1);
        //подсчёт значений для гистограммы
        for(int i = 0;i<val.size();i++) {
            if(val[i] > ranges[k+1]) {
                nums[k]/=(val.size()*(ranges[k+1]-ranges[k]));
                while(val[i]>ranges[k+1])
                    k++;
            }
            nums[k]++;
        }
        nums[nums.size()-1]/=(val.size()*ranges[k+1]-ranges[k]);
        valX.resize(3*(ranges.size()-1)+1); //для красивой гистограммы
        hist.resize(3*(ranges.size()-1)+1); //
        ui->histStats->setColumnCount(ranges.size()-1);

        //заполнение таблицы
        for(int i = 0;i<ranges.size()-1;i++){
            double midx = (ranges[i+1]+ranges[i])/2;

            valX[3*i] = valX[3*i+1] = ranges[i];
            valX[3*i+2] = ranges[i+1];
            hist[3*i] = 0;
            hist[3*i+1] = hist[3*i+2] = nums[i];

            tbl = new QTableWidgetItem(QString::number(midx));
            ui->histStats->setItem(0,i,tbl);
            tbl = new QTableWidgetItem(QString::number(densf(lambda,midx)));
            ui->histStats->setItem(1,i,tbl);
            tbl = new QTableWidgetItem(QString::number(nums[i]));
            ui->histStats->setItem(2,i,tbl);
            if(maxdif<std::abs(nums[i]-densf(lambda,midx)))
                maxdif = std::abs(nums[i]-densf(lambda,midx));
        }
        valX[valX.size()-1] = ranges[ranges.size()-1];
        hist[hist.size()-1] = 0;

        paintGraph(ui->histGraph,unX,f,"f(x)",valX,hist,"Гистогр.");
        ui->max->setText(QString::number(maxdif));
    }

}


//При изменении количества отрезков (или выборки)
void MainWindow::on_histN_textChanged(const QString &arg)
{
    QTableWidgetItem* tbl;
    int num = arg.toInt();
    ui->histRanges->setColumnCount(num+1);
    if(!val.isEmpty()){

        tbl = new QTableWidgetItem(QString::number(val[0]));
        ui->histRanges->setItem(0,0,tbl);
        tbl = new QTableWidgetItem(QString::number(val[val.size()-1]));
        ui->histRanges->setItem(0,num,tbl);
        for(int i = 1; i <num;i++){
            tbl = new QTableWidgetItem(QString::number(val[0]+i*(val[val.size()-1]-val[0])/num));
            ui->histRanges->setItem(0,i,tbl);

        }
    }
}


//При изменении кол-ва отрезков (или выборки)
void MainWindow::on_hypN_textChanged(const QString &arg)
{
    ui->hypRanges->blockSignals(true);
    QTableWidgetItem* tbl;
    int num = arg.toInt();
    ui->hypRanges->setColumnCount(num+1);
    ui->hypH0->setColumnCount(num);
    double lambda = ui->lambdaPar->text().toDouble();

    tbl = new QTableWidgetItem(QString::number(-std::numeric_limits<double>::infinity()));
    ui->hypRanges->setItem(0,0,tbl);
    tbl = new QTableWidgetItem(QString::number(std::numeric_limits<double>::infinity()));
    ui->hypRanges->setItem(0,num,tbl);

    double beg,end;

    if(!val.isEmpty()){
        beg = val[0];
        end = val[val.size()-1];
    } else {
        beg = 0;
        end = 0;
    }

    double r,ro = -std::numeric_limits<double>::infinity();
    for(int i = 1; i <num;i++){
        r = beg+i*(end-beg)/num;

        tbl = new QTableWidgetItem(QString::number(intF(lambda,r) - intF(lambda,ro)));
        ui->hypH0->setItem(0,i-1,tbl);

        ro = r;

        tbl = new QTableWidgetItem(QString::number(r));
        ui->hypRanges->setItem(0,i,tbl);
    }
    tbl = new QTableWidgetItem(QString::number(1-intF(lambda,r)));
    ui->hypH0->setItem(0,num-1,tbl);
    ui->hypRanges->blockSignals(false);
}

double chi(double r,double x) {
    if(x > 0)
        return std::pow(x,r/2-1)*exp(-x/2)/tgamma(r/2)/std::pow(2,r/2);
    else
        return 0;

}

double Frev(double r, double x){
    int n = (int)(x*1000);
    double h = x/n;
    double res=0;
    for(int i = 0; i < n; i++)
        res+=(chi(r,i*h)+chi(r,(i+1)*h))*h/2;
    return 1-res;
}

//Проверка гипотезы
void MainWindow::on_hypStart_clicked()
{
    if(!val.isEmpty()){

        QVector<double> q,ranges,nums;
        //Вытащить границы разбиения из таблицы
        ranges.push_back(-std::numeric_limits<double>::infinity());
        q.push_back(ui->hypH0->item(0,0)->text().toDouble());

        for(int i = 1; i < ui->hypRanges->columnCount()-1; i++) {
            ranges.push_back(ui->hypRanges->item(0,i)->text().toDouble());
            q.push_back(ui->hypH0->item(0,i)->text().toDouble());
        }
        ranges.push_back(std::numeric_limits<double>::infinity());

        int k = 0;
        nums.resize(ranges.size()-1);

        //подсчёт значений для R0

        for(int i = 0;i<val.size();i++) {
            if(val[i] > ranges[k+1]) {
                while(val[i]>ranges[k+1])
                    k++;
            }
            nums[k]++;
        }

        //Подсчёт R0 и Frev(R0)
        double r0 = 0;
        int n = ui->nPar->text().toInt();
        for(int i = 0; i < q.size(); i++)
            r0 += (nums[i] - n*q[i])*(nums[i] - n*q[i])/(n*q[i]);


        double testval = Frev(q.size()-1,r0);

        //Вывод
        if(testval>ui->hypSign->text().toDouble())
            ui->hypRes->setText("Гипотеза принята");
        else
            ui->hypRes->setText("Гипотеза отвергнута");
        ui->hypTest->setText(QString::number(testval));
    }
}



void MainWindow::on_hypRanges_itemChanged(QTableWidgetItem *item)
{

    QTableWidgetItem* tbl;
    double lambda = ui->lambdaPar->text().toDouble();
    for(int i = 0; i < ui->hypH0->columnCount();i++){
        tbl = new QTableWidgetItem(QString::number(intF(lambda,ui->hypRanges->item(0,i+1)->text().toDouble()) - intF(lambda,ui->hypRanges->item(0,i)->text().toDouble())));
        ui->hypH0->setItem(0,i,tbl);
    }
}
