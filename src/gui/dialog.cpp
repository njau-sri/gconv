#include <QDir>
#include <QDateTime>
#include <QFileDialog>
#include <QMessageBox>
#include "dialog.h"
#include "ui_dialog.h"

Dialog::Dialog(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::Dialog),
    proc_(0)
{
    ui->setupUi(this);

    proc_ = new QProcess(this);
    proc_->setProcessChannelMode(QProcess::MergedChannels);
    connect(proc_, SIGNAL(readyReadStandardOutput()), this, SLOT(slot_proc_readyReadStandardOutput()));
    connect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));

    setWindowTitle(QLatin1String("GCONV 2019.0.dev"));
}

Dialog::~Dialog()
{
    delete ui;
}

void Dialog::closeEvent(QCloseEvent *e)
{
    if ( isProcessRunning() ) {
        e->ignore();
        return;
    }

    e->accept();
}

void Dialog::on_pushButtonIn_clicked()
{
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose Input Genotype File"));
    if ( ! fileName.isEmpty() )
        ui->lineEditIn->setText(fileName);
}

void Dialog::on_pushButtonOut_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Choose Output Genotype File"));
    if ( ! fileName.isEmpty() )
        ui->lineEditOut->setText(fileName);
}

void Dialog::on_buttonBox_accepted()
{
    QString prog = QDir(QApplication::applicationDirPath()).filePath(QLatin1String("gconv"));

    QStringList format, suffix;

    format << QLatin1String("--vcf")
           << QLatin1String("--ped")
           << QLatin1String("--hmp")
           << QLatin1String("--geno");

    suffix << QLatin1String(".vcf")
           << QLatin1String(".ped")
           << QLatin1String(".hmp")
           << QLatin1String(".geno");

    QStringList args;

    if ( ! ui->lineEditIn->text().isEmpty() )
        args << format.at(ui->comboBoxIn->currentIndex()) << ui->lineEditIn->text();

    if ( ! ui->lineEditOut->text().isEmpty() )
        args << QLatin1String("--out") << ui->lineEditOut->text() + suffix.at(ui->comboBoxOut->currentIndex());

    if ( ui->checkBoxSort->isChecked() )
        args << QLatin1String("--sort");

    startProcess(prog, args);
}

void Dialog::on_buttonBox_rejected()
{
    if (ui->buttonBox->button(QDialogButtonBox::Ok)->isEnabled())
        close();

    isProcessRunning();
}

void Dialog::slot_proc_readyReadStandardOutput()
{
    QByteArray v = proc_->readAllStandardOutput();
    ui->plainTextEditLog->moveCursor(QTextCursor::End);
    ui->plainTextEditLog->insertPlainText(QString::fromLocal8Bit(v.data(),v.size()));
    ui->plainTextEditLog->moveCursor(QTextCursor::End);
}

void Dialog::slot_proc_finished(int code, QProcess::ExitStatus status)
{
    ui->progressBar->setRange(0,1);
    ui->progressBar->reset();
    ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);

    if (status != QProcess::NormalExit || code != 0) {
        QMessageBox::critical(this, tr("ERROR"), tr("Process exited unexpectedly: code %1, status %2.").arg(code).arg(status));
        return;
    }
}

bool Dialog::isProcessRunning()
{
    if (proc_->state() != QProcess::NotRunning) {
        if (QMessageBox::question(this, tr("Terminate"), tr("Stop currently running computations?"), QMessageBox::Ok | QMessageBox::Cancel) != QMessageBox::Ok)
            return true;

        disconnect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));
        proc_->close();

        ui->progressBar->setRange(0,1);
        ui->progressBar->reset();
        ui->buttonBox->button(QDialogButtonBox::Ok)->setEnabled(true);

        connect(proc_, SIGNAL(finished(int,QProcess::ExitStatus)), this, SLOT(slot_proc_finished(int,QProcess::ExitStatus)));
    }

    return false;
}

void Dialog::startProcess(const QString &prog, const QStringList &args)
{
    if ( isProcessRunning() )
        return;

    proc_->start(prog, args);
    if ( ! proc_->waitForStarted() ) {
        QMessageBox::critical(this, tr("ERROR"), tr("Can't start process: %1").arg(prog));
        return;
    }

    ui->progressBar->setRange(0,0);
    ui->progressBar->reset();
    ui->buttonBox->button(QDialogButtonBox::Ok)->setDisabled(true);
}
