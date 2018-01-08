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
    QString fileName = QFileDialog::getOpenFileName(this, tr("Choose input genotype file"));
    if (!fileName.isEmpty())
        ui->lineEditIn->setText(fileName);
}

void Dialog::on_pushButtonOut_clicked()
{
    QString fileName = QFileDialog::getSaveFileName(this, tr("Choose output genotype file"));
    if (!fileName.isEmpty())
        ui->lineEditOut->setText(fileName);
}

void Dialog::on_buttonBox_accepted()
{
    QString prog = QDir(QApplication::applicationDirPath()).filePath(QLatin1String("gconv"));

    QStringList args;

    if (ui->comboBoxIn->currentIndex() != 0 && ! ui->lineEditIn->text().isEmpty()) {
        if (ui->comboBoxIn->currentIndex() == 1)
            args << QLatin1String("--vcf") << ui->lineEditIn->text();
        else if (ui->comboBoxIn->currentIndex() == 2)
            args << QLatin1String("--ped") << ui->lineEditIn->text();
        else if (ui->comboBoxIn->currentIndex() == 3)
            args << QLatin1String("--hmp") << ui->lineEditIn->text();
        else if (ui->comboBoxIn->currentIndex() == 4)
            args << QLatin1String("--geno") << ui->lineEditIn->text();
    }

    if (ui->comboBoxOut->currentIndex() != 0 && ! ui->lineEditOut->text().isEmpty()) {
        QString suffix;
        if (ui->comboBoxOut->currentIndex() == 1)
            suffix = QLatin1String(".vcf");
        else if (ui->comboBoxOut->currentIndex() == 2)
            suffix = QLatin1String(".ped");
        else if (ui->comboBoxOut->currentIndex() == 3)
            suffix = QLatin1String(".hmp");
        else if (ui->comboBoxOut->currentIndex() == 4)
            suffix = QLatin1String(".geno");
        QString out = ui->lineEditOut->text();
        if ( ! out.endsWith(suffix) )
            out.append(suffix);
        args << QLatin1String("--out") << out;
    }

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
