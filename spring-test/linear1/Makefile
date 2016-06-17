
# Fortran コンパイラの指定

FC=ifort

# Fortran オプション

FFLAGS=/nologo
FFLAGS=$(FFLAGS) /recursive     # compile all procedures for possible recursive execution (再帰可能関数としてコンパイル)
FFLAGS=$(FFLAGS) /Qauto-scalar  # make scalar local variables AUTOMATIC  (saveにしない）
FFLAGS=$(FFLAGS) /QaxWPT        # Core2Duo, Core2Quad, Xeon用でSSSE3に対応したコードを作成する．
FFLAGS=$(FFLAGS) /traceback     # 実行時エラー時で，数値でなくシンボル名でトレースバックを表示する為のPC correlationデータを作成する．
FFLAGS=$(FFLAGS) /O2            # 速度最適化
FFLAGS=$(FFLAGS) /Qzero         # 初期値が設定されていないsave局所変数のゼロ初期化を行う
FFLAGS=$(FFLAGS) /Qautodouble   # set default size of REAL to 8 bytes
FFLAGS=$(FFLAGS) /RTCu          # 初期化されていない変数の使用をレポートする
FFLAGS=$(FFLAGS) /Zi            # オブジェクトファイルにデバッグシンボルを作成．


# ライブラリ作成オプション

LIBFLAGS=/nologo
LIBFLAGS=$(LIBFLAGS) /subsystem:console  # コンソールプログラム
LIBFLAGS=$(LIBFLAGS) /machine:X64    # 64bit machine用
# LIBFLAGS=$(LIBFLAGS) /NODEFAULTLIB:libifcoremt.lib  

# リンクフラグ

LINKFLAGS=$(LIBFLAGS)
LINKFLAGS=$(LINKFLAGS) /INCREMENTAL:NO
#LINKFLAGS=$(LINKFLAGS) /NODEFAULTLIB:LIBC.LIB
#LINKFLAGS=$(LINKFLAGS) /NODEFAULTLIB:LIBCMT.LIB
#LINKFLAGS=$(LINKFLAGS) /DEFAULTLIB:OLDNAMES.LIB
#LINKFLAGS=$(LINKFLAGS) /DEFAULTLIB:LIBIFCOREMD.LIB
#LINKFLAGS=$(LINKFLAGS) /DEFAULTLIB:LIBIFPORTMD.LIB
#LINKFLAGS=$(LINKFLAGS) /DEFAULTLIB:LIBMMD.LIB
#LINKFLAGS=$(LINKFLAGS) /DEFAULTLIB:MSVCRT.LIB
#LINKFLAGS=$(LINKFLAGS) /DEFAULTLIB:kernel32.lib
#LINKFLAGS=$(LINKFLAGS) /DEFAULTLIB:user32.lib
#LINKFLAGS=$(LINKFLAGS) /DEFAULTLIB:advapi32.lib
#LINKFLAGS=$(LINKFLAGS) $(MKL)
#LINKFLAGS=$(LINKFLAGS) $(CXX_LIB)
LINKFLAGS=$(LINKFLAGS) /FIXED:NO
LINKFLAGS=$(LINKFLAGS) /dll


# DLLからエクスポートするシンボルの指定．(これはKから始まっている必要がある）
EXPORTS=/EXPORT:K_LINEAR_SPRING


o=obj
FMODULE=debug.$o linear_spring.$o bilinear_spring.$o
FOBJS=extdb.$o

run:check.dat
	echo 1 | dat_extract check.dat


check.dat: exec.mak abaqus_v6.env check.inp
	nmake -f exec.mak

exec.mak:mkmake.exe
	dir /b *.inp | mkmake.exe int > exec.mak

abaqus_v6.env: env_create.exe uels.lib
	env_create.exe

all:lib

lib:uels.lib

uels.lib:$(FOBJS) $(FMODULE)
	LIB $(LIBFLAGS) /OUT:$@ $(FOBJS) $(FMODULE)

clean:
	del $(FOBJS) $(FMODULE) *.mod

.f90.$o:
	$(FC) /c $(FFLAGS) $<


.f.$o:
	$(FC) /c $(FFLAGS) $<


nonlinear_spring.$o:debug.$o
extdb.$o:debug.$o

.PHONY: clean lib run

