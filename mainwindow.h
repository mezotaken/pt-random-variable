#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <qcustomplot.h>
#include <ctime>
#include <cmath>
#include <random>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_startButton_clicked();
    void paintGraph(QCustomPlot* dst,QVector<double>&X,QVector<double>& F1,QString name1, QVector<double>& Xaux,QVector<double>& F2,QString name2);

    void on_isrdSeed_stateChanged(int arg1);

    void on_histStart_clicked();

    void on_histN_textChanged(const QString &num);

    void on_hypN_textChanged(const QString &arg1);

    void on_hypStart_clicked();

    void on_hypRanges_itemChanged(QTableWidgetItem *item);

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
