import pytest


@pytest.fixture(scope='module', params=["classifier", "regressor"])
def boosting_alg_type(request):
    return request.param


@pytest.fixture(scope='module',
    params=["AdaBoost", "GBoost", "LR", "NGBoost", "RF", "SGBoost", "SVM_RBF", "SVML", "XGBoost"])
def boosting_method(request):
    return request.param
