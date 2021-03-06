#ifndef DIALOG_H
#define DIALOG_H

#include <QDialog>
#include <QProcess>

namespace Ui {
class Dialog;
}

class Dialog : public QDialog
{
    Q_OBJECT

public:

    explicit Dialog(QWidget *parent = 0);

    ~Dialog();

protected:

    void closeEvent(QCloseEvent *e);

private slots:

    void on_pushButtonIn_clicked();

    void on_pushButtonOut_clicked();

    void on_buttonBox_accepted();

    void on_buttonBox_rejected();

    void slot_proc_readyReadStandardOutput();

    void slot_proc_finished(int code, QProcess::ExitStatus status);

private:

    bool isProcessRunning();

    void startProcess(const QString &prog, const QStringList &args);

private:

    Ui::Dialog *ui;
    QProcess *proc_;
};

#endif // DIALOG_H
